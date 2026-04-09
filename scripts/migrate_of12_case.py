#!/usr/bin/env python3
"""
OpenFOAM-8 to OpenFOAM-12 case migration helper.

Default mode is dry-run. Pass --apply to write changes.

This script is intentionally conservative. It applies the mechanical changes
used for Mesh_3x3 and reports remaining items that need manual review.
"""

from __future__ import annotations

import argparse
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path


RAW_TABLE_FILES = {"IDN", "Idif", "sunPosVector"}
BACKUP_SUFFIXES = (".bak", ".bak_of8")
NON_REGION_SYSTEM_FILES = {
    "blockMeshDict",
    "blockMeshDict.air",
    "blockMeshDict.veg",
    "controlDict",
    "createPatchDict",
    "decomposeParDict",
    "fvSchemes",
    "fvSolution",
    "meshQualityDict",
    "pedestrianCloudPoints",
    "snappyHexMeshDict",
    "snappyHexMeshDict.air",
    "snappyHexMeshDict.veg",
    "surfaceFeaturesDict",
    "viewFactorsDict",
}


@dataclass
class Report:
    changed: list[str] = field(default_factory=list)
    review: list[str] = field(default_factory=list)
    unchanged: list[str] = field(default_factory=list)


@dataclass
class RegionInfo:
    name: str
    kind: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Migrate an urbanMicroclimateFoam case to OpenFOAM-12 syntax."
    )
    parser.add_argument(
        "case",
        type=Path,
        nargs="?",
        default=Path("/home/strebdom/Mesh_3x3"),
        help="Case directory to migrate. Default: /home/strebdom/Mesh_3x3",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write changes. Without this flag, only print what would change.",
    )
    parser.add_argument(
        "--keep-fvoptions",
        action="store_true",
        help="Do not rename system/<region>/fvOptions after migrating constraints.",
    )
    return parser.parse_args()


def is_backup(path: Path) -> bool:
    return path.name.endswith(BACKUP_SUFFIXES) or ".bak_" in path.name


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="surrogateescape")


def write_text(path: Path, text: str, apply: bool) -> None:
    if apply:
        path.write_text(text, encoding="utf-8")


def update_file(path: Path, transform, report: Report, apply: bool) -> None:
    original = read_text(path)
    updated = transform(original)
    if updated != original:
        write_text(path, updated, apply)
        report.changed.append(str(path))


def all_files(case_dir: Path):
    for root in ("system", "constant", "0"):
        base = case_dir / root
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and not is_backup(path):
                yield path


def parse_region_groups(text: str) -> list[tuple[str, list[str]]]:
    region_groups: list[tuple[str, list[str]]] = []
    for key in ("fluid", "solid", "vegetation"):
        match = re.search(rf"{key}\s*\(([^)]*)\)", text, flags=re.S)
        if match and match.group(1).strip():
            region_groups.append((key, match.group(1).split()))
    return region_groups


def detect_regions(case_dir: Path) -> list[RegionInfo]:
    regions: dict[str, RegionInfo] = {}

    props_path = case_dir / "constant/regionProperties"
    if props_path.exists():
        for kind, names in parse_region_groups(read_text(props_path)):
            for name in names:
                regions[name] = RegionInfo(name=name, kind=kind)

    system_dir = case_dir / "system"
    if system_dir.exists():
        for path in sorted(system_dir.iterdir()):
            if not path.is_dir() or path.name in NON_REGION_SYSTEM_FILES:
                continue
            info = regions.get(path.name, RegionInfo(name=path.name))
            regions[path.name] = info

    for root in ("constant", "0"):
        base = case_dir / root
        if not base.exists():
            continue
        for path in sorted(base.iterdir()):
            if not path.is_dir():
                continue
            if root == "constant" and path.name in {"extendedFeatureEdgeMesh", "triSurface"}:
                continue
            info = regions.get(path.name, RegionInfo(name=path.name))
            regions[path.name] = info

    return [regions[name] for name in sorted(regions)]


def parse_number_of_subdomains(path: Path) -> int | None:
    if not path.exists():
        return None
    match = re.search(r"\bnumberOfSubdomains\s+(\d+)\s*;", read_text(path))
    return int(match.group(1)) if match else None


def shell_quote_words(words: list[str]) -> str:
    return " ".join(words)


def run_application_line(cmd: str, args: str = "") -> str:
    suffix = f" {args}" if args else ""
    return f'runApplication {cmd}{suffix}'


def logged_command_line(cmd: str, log_name: str, args: str = "") -> str:
    suffix = f" {args}" if args else ""
    return f"{cmd}{suffix} > {log_name} 2>&1"


def stage_region_poly_mesh_lines(region: str) -> list[str]:
    return [
        f"if [ -d constant/polyMesh ]; then",
        f"    rm -rf constant/{region}/polyMesh",
        f"    mv constant/polyMesh constant/{region}/.",
        "fi",
        f'if [ ! -d constant/{region}/polyMesh ]; then',
        f'    echo "Missing polyMesh for region {region}" >&2',
        "    exit 1",
        "fi",
    ]


def region_mesh_variant(case_dir: Path, region: str, kind: str | None) -> tuple[str, str] | None:
    candidates: list[tuple[str, str]] = [
        (f"blockMeshDict.{region}", f"snappyHexMeshDict.{region}"),
    ]
    if kind == "fluid":
        candidates.append(("blockMeshDict.air", "snappyHexMeshDict.air"))
    elif kind == "vegetation":
        candidates.append(("blockMeshDict.veg", "snappyHexMeshDict.veg"))

    for block_name, snappy_name in candidates:
        if (case_dir / "system" / block_name).exists() and (case_dir / "system" / snappy_name).exists():
            return block_name, snappy_name
    return None


def generated_allprepare(case_dir: Path) -> str:
    regions = detect_regions(case_dir)
    fluid_regions = [r.name for r in regions if r.kind == "fluid"]
    vegetation_regions = [r.name for r in regions if r.kind == "vegetation"]
    solid_regions = [r.name for r in regions if r.kind == "solid"]
    other_regions = [r.name for r in regions if r.kind not in {"fluid", "solid", "vegetation"}]

    lines = [
        "#!/bin/bash",
        "cd ${0%/*} || exit 1",
        "",
        '. "$WM_PROJECT_DIR/bin/tools/RunFunctions"',
        "",
        "set -e",
        "",
        "rm -rf log",
        "",
    ]

    if (case_dir / "system/surfaceFeaturesDict").exists():
        lines.append(run_application_line("surfaceFeatures"))
        lines.append("")

    staged_regions = fluid_regions + vegetation_regions
    for region in staged_regions:
        info = next(r for r in regions if r.name == region)
        variant = region_mesh_variant(case_dir, region, info.kind)
        if not variant:
            continue

        block_name, snappy_name = variant
        lines.extend(
            [
                f'echo "Creating mesh for {region} region"',
                f"cp system/{block_name} system/blockMeshDict",
                logged_command_line("blockMesh", f"log.blockMesh.{region}"),
            ]
        )
        if (case_dir / "system/createPatchDict").exists():
            lines.append(logged_command_line("createPatch", f"log.createPatch.{region}", "-overwrite"))
        lines.extend(
            [
                f"cp system/{snappy_name} system/snappyHexMeshDict",
                logged_command_line("snappyHexMesh", f"log.snappyHexMesh.{region}", "-overwrite"),
            ]
        )
        lines.extend(stage_region_poly_mesh_lines(region))
        if (case_dir / f"system/{region}/changeDictionaryDict").exists():
            lines.append(f"changeDictionary -region {region} > log.changeDictionary.{region} 2>&1")
        lines.append("")

        if (case_dir / f"system/{region}/topoSetDict").exists():
            lines.extend(
                [
                    f'echo "Running topoSet for {region} region"',
                    f"topoSet -region {region} > log.topoSet.{region} 2>&1",
                    "",
                ]
            )

        if (case_dir / f"system/{region}/setFieldsDict").exists():
            lad_backup = case_dir / f"0/{region}.bckp/LAD"
            lines.extend([f'echo "Setting fields for {region} region"'])
            if lad_backup.exists():
                lines.extend(
                    [
                        f"rm -f 0/{region}/LAD 0/{region}/LAD.gz",
                        f"cp 0/{region}.bckp/LAD 0/{region}/.",
                    ]
                )
            lines.extend(
                [
                    f"setFields -region {region} > log.setFields.{region} 2>&1",
                    "",
                ]
            )

    loop_regions = [
        name
        for name in solid_regions + other_regions
        if (case_dir / f"system/{name}/blockMeshDict").exists()
    ]
    if loop_regions:
        lines.extend(
            [
                f'echo "Creating mesh for {" ".join(loop_regions)} regions"',
                f"for i in {shell_quote_words(loop_regions)}",
                "do",
                "    blockMesh -region $i > log.blockMesh.$i 2>&1",
            ]
        )
        if any((case_dir / f"system/{name}/createPatchDict").exists() for name in loop_regions):
            lines.append("    if [ -f system/$i/createPatchDict ]; then createPatch -region $i -overwrite; fi")
        if any((case_dir / f"system/{name}/changeDictionaryDict").exists() for name in loop_regions):
            lines.append("    if [ -f system/$i/changeDictionaryDict ]; then changeDictionary -region $i > log.changeDictionary.$i 2>&1; fi")
        if any((case_dir / f"system/{name}/topoSetDict").exists() for name in loop_regions):
            lines.append("    if [ -f system/$i/topoSetDict ]; then topoSet -region $i > log.topoSet.$i 2>&1; fi")
        lines.extend(["done", ""])

    lines.extend(
        [
            "mkdir -p log",
            'find . -maxdepth 1 -type f -name "log.*" -exec mv {} log/. \\;',
            "",
            "# ----------------------------------------------------------------- end-of-file",
            "",
        ]
    )
    return "\n".join(lines)


def generated_allrun(case_dir: Path) -> str:
    regions = detect_regions(case_dir)
    fluid_regions = [r.name for r in regions if r.kind == "fluid"]
    vegetation_regions = [r.name for r in regions if r.kind == "vegetation"]
    nprocs = (
        parse_number_of_subdomains(case_dir / "system/decomposeParDict")
        or next(
            (
                parse_number_of_subdomains(case_dir / f"system/{region}/decomposeParDict")
                for region in [r.name for r in regions]
                if parse_number_of_subdomains(case_dir / f"system/{region}/decomposeParDict") is not None
            ),
            None,
        )
        or 20
    )

    lines = [
        "#!/bin/bash",
        "cd ${0%/*} || exit 1",
        "",
        '. "$WM_PROJECT_DIR/bin/tools/RunFunctions"',
        "",
        "set -e",
        "",
        "rm -rf processor*",
        "",
    ]

    for region in vegetation_regions:
        source = case_dir / f"constant/{region}/viewFactorsDict"
        if source.exists():
            lines.extend(
                [
                    f"cp constant/{region}/viewFactorsDict system/viewFactorsDict",
                    f"cp constant/{region}/viewFactorsDict system/{region}/viewFactorsDict",
                    "",
                ]
            )

    lines.append("decomposePar -allRegions")
    lines.append("")

    staged_region_files: list[tuple[str, str]] = []
    for region in fluid_regions + vegetation_regions:
        for name in ("radiationProperties", "solarLoadProperties"):
            if (case_dir / f"constant/{region}/{name}").exists():
                staged_region_files.append((region, name))
    for region in fluid_regions:
        if (case_dir / f"constant/{region}/vegetationProperties").exists():
            staged_region_files.append((region, "vegetationProperties"))

    if staged_region_files:
        lines.extend(
            [
                "for p in processor*; do",
                '    [ -d "$p" ] || continue',
            ]
        )
        for region, name in staged_region_files:
            lines.append(f'    mkdir -p "$p/constant/{region}"')
            lines.append(f'    cp constant/{region}/{name} "$p/constant/{region}/"')
        lines.extend(["done", ""])

    for region in vegetation_regions:
        if (case_dir / f"constant/{region}/viewFactorsDict").exists():
            lines.append(f"mpirun -np {nprocs} faceAgglomerate -region {region} -parallel")

    for region in fluid_regions:
        if (case_dir / f"0/{region}/LAD").exists() or (case_dir / f"0/{region}/LAD.gz").exists():
            lines.append(f"mpirun -np {nprocs} calcLAI -region {region} -parallel")

    for region in vegetation_regions:
        if (case_dir / f"constant/{region}/viewFactorsDict").exists():
            lines.append(f"mpirun -np {nprocs} viewFactorsGen -region {region} -parallel")
            lines.append(f"mpirun -np {nprocs} solarRayTracingGen -region {region} -parallel")

    lines.extend(
        [
            'echo "Running solver"',
            f"mpirun -np {nprocs} urbanMicroclimateFoam -parallel > log.UMC 2>&1",
            "reconstructPar -allRegions",
            "",
            "# ----------------------------------------------------------------- end-of-file",
            "",
        ]
    )
    return "\n".join(lines)


def generated_allclean(case_dir: Path) -> str:
    regions = detect_regions(case_dir)
    lines = [
        "#!/bin/bash",
        "cd ${0%/*} || exit 1",
        "",
        '. "$WM_PROJECT_DIR/bin/tools/CleanFunctions"',
        "",
        "cleanCase",
        "",
    ]

    for region in [r.name for r in regions]:
        if (case_dir / "constant" / region).exists():
            lines.append(f"foamCleanPolyMesh -region {region}")

    radiation_regions = [
        r.name
        for r in regions
        if (case_dir / f"constant/{r.name}/viewFactorsDict").exists()
    ]
    if radiation_regions:
        lines.append("")
        for region in radiation_regions:
            lines.extend(
                [
                    f"rm -f constant/{region}/F*",
                    f"rm -f constant/{region}/constructMap*",
                    f"rm -f constant/{region}/finalAgglom*",
                    f"rm -f constant/{region}/globalFaceFaces*",
                    f"rm -f constant/{region}/subMap*",
                    f"rm -f constant/{region}/skyViewCoeff*",
                    f"rm -f constant/{region}/sunViewCoeff*",
                    f"rm -f constant/{region}/sunskyMap*",
                    f"rm -f constant/{region}/sunVisibleOrNot*",
                ]
            )
            if (case_dir / "0" / region).exists():
                lines.append(f"rm -f 0/{region}/viewFactorField*")

    lines.extend(
        [
            "",
            "# ----------------------------------------------------------------- end-of-file",
            "",
        ]
    )
    return "\n".join(lines)


def ensure_executable(path: Path, apply: bool) -> None:
    if apply and path.exists():
        path.chmod(path.stat().st_mode | 0o111)


def ensure_run_scripts(case_dir: Path, report: Report, apply: bool) -> None:
    scripts = {
        case_dir / "Allprepare": generated_allprepare(case_dir),
        case_dir / "Allrun": generated_allrun(case_dir),
        case_dir / "Allclean": generated_allclean(case_dir),
    }

    for path, updated in scripts.items():
        original = read_text(path) if path.exists() else None
        if original != updated:
            write_text(path, updated, apply)
            ensure_executable(path, apply)
            report.changed.append(str(path))


def migrate_momentum_transport(text: str) -> str:
    text = re.sub(r"object\s+RASProperties;", "object      momentumTransport;", text)
    text = re.sub(r"(^\s*)RASModel(\s+\S+\s*;)", r"\1model\2", text, flags=re.M)
    return text


def migrate_fvsolution(text: str) -> str:
    text = re.sub(r"^\s*pMaxFactor\s+[^;]+;\s*\n", "", text, flags=re.M)
    text = re.sub(r"^\s*pMinFactor\s+[^;]+;\s*\n", "", text, flags=re.M)
    return text


def migrate_setfields(text: str) -> str:
    lines = text.splitlines(keepends=True)
    out: list[str] = []
    in_zone_to_cell = False
    depth = 0
    pending_zone_to_cell = False

    for line in lines:
        if re.match(r"^\s*zoneToCell\s*$", line):
            pending_zone_to_cell = True

        if pending_zone_to_cell and "{" in line:
            in_zone_to_cell = True
            pending_zone_to_cell = False

        if in_zone_to_cell:
            line = re.sub(r"(^\s*)name(\s+\S+\s*;)", r"\1zone\2", line)

        out.append(line)

        if in_zone_to_cell or pending_zone_to_cell:
            depth += line.count("{") - line.count("}")
            if in_zone_to_cell and depth <= 0:
                in_zone_to_cell = False
                depth = 0

    return "".join(out)


def migrate_blockmesh(text: str) -> str:
    return text.replace("$:", "$!")


def transform_patch_block(block: str) -> str:
    had_nearest_cell = re.search(r"sampleMode\s+nearestCell\s*;", block) is not None

    block = re.sub(r"(^\s*)type(\s+)mappedPatch\s*;", r"\1type\2mapped;", block, flags=re.M)
    block = re.sub(r"(^\s*)sampleRegion(\s+)", r"\1neighbourRegion\2", block, flags=re.M)
    block = re.sub(r"(^\s*)samplePatch(\s+)", r"\1neighbourPatch\2", block, flags=re.M)
    block = re.sub(r"(^\s*)sampleMode\s+nearest\s*;", r"\1method          nearest;", block, flags=re.M)
    block = re.sub(r"^\s*sampleMode\s+nearestPatchFace\s*;\s*\n", "", block, flags=re.M)

    if had_nearest_cell:
        block = re.sub(r"^\s*sampleMode\s+nearestCell\s*;\s*\n", "", block, flags=re.M)
        block = re.sub(r"(^\s*)type\s+(mappedWall|mappedPatch|mapped)\s*;", r"\1type            mappedInternal;", block, flags=re.M)
        block = re.sub(r"(^\s*)offsetMode\s+\S+\s*;", r"\1offsetMode      normal;", block, flags=re.M)
        block = re.sub(r"^\s*offset\s+\([^;]+\)\s*;\s*\n", "", block, flags=re.M)
        if "distance" not in block:
            block = re.sub(r"(^\s*offsetMode\s+normal\s*;\s*\n)", r"\1        distance        0;\n", block, flags=re.M)

    return block


def migrate_change_dictionary(text: str) -> str:
    lines = text.splitlines(keepends=True)
    out: list[str] = []
    block: list[str] = []
    in_patch = False
    pending_patch_name = False
    depth = 0

    for line in lines:
        if not in_patch:
            if re.match(r"^\s*[A-Za-z0-9_./:-]+\s*$", line):
                pending_patch_name = True
                block = [line]
                continue

            if pending_patch_name and "{" in line:
                in_patch = True
                depth = line.count("{") - line.count("}")
                block.append(line)
                if depth <= 0:
                    out.append(transform_patch_block("".join(block)))
                    block = []
                    in_patch = False
                    pending_patch_name = False
                continue

            if pending_patch_name:
                out.extend(block)
                block = []
                pending_patch_name = False

            out.append(line)
            continue

        block.append(line)
        depth += line.count("{") - line.count("}")
        if depth <= 0:
            out.append(transform_patch_block("".join(block)))
            block = []
            in_patch = False
            pending_patch_name = False

    if block:
        out.extend(block)

    return "".join(out)


def mapped_internal_patches(boundary_text: str) -> set[str]:
    matches = re.finditer(
        r"^\s*([A-Za-z0-9_./:-]+)\s*\n\s*\{\s*\n(?:.*\n)*?\s*type\s+mappedInternal\s*;",
        boundary_text,
        flags=re.M,
    )
    return {match.group(1) for match in matches}


def migrate_mapped_internal_field(text: str, patch_names: set[str]) -> str:
    if not patch_names or "mapped;" not in text:
        return text

    updated = text
    for patch_name in sorted(patch_names):
        updated = re.sub(
            rf"(^\s*{re.escape(patch_name)}\s*\n\s*\{{.*?^\s*type\s+)mapped(\s*;)",
            rf"\1mappedInternalValue\2",
            updated,
            flags=re.M | re.S,
        )

    return updated


def migrate_region_initial_fields(case_dir: Path, report: Report, apply: bool) -> None:
    for region_info in detect_regions(case_dir):
        region = region_info.name
        field_dir = case_dir / "0" / region
        if not field_dir.is_dir():
            continue

        patches: set[str] = set()

        boundary_path = case_dir / "constant" / region / "polyMesh/boundary"
        if boundary_path.exists():
            patches = mapped_internal_patches(read_text(boundary_path))

        if not patches:
            change_dict_path = case_dir / "system" / region / "changeDictionaryDict"
            if change_dict_path.exists():
                patches = mapped_internal_patches(read_text(change_dict_path))

        if not patches:
            continue

        for field_path in sorted(field_dir.iterdir()):
            if not field_path.is_file() or is_backup(field_path):
                continue
            update_file(
                field_path,
                lambda text, names=patches: migrate_mapped_internal_field(text, names),
                report,
                apply,
            )


def regions_and_region_solvers_from_region_properties(path: Path) -> tuple[str, str] | None:
    if not path.exists():
        return None

    text = read_text(path)
    region_groups: list[tuple[str, list[str]]] = []
    solver_entries: list[tuple[str, str]] = []
    for key in ("fluid", "solid", "vegetation"):
        match = re.search(rf"{key}\s*\(([^)]*)\)", text, flags=re.S)
        if match and match.group(1).strip():
            names = match.group(1).split()
            region_groups.append((key, names))
            for name in names:
                solver_entries.append((name, key))

    if not region_groups:
        return None

    regions_lines = ["regions", "{"]
    for key, names in region_groups:
        joined = " ".join(names)
        regions_lines.append(f"    {key:<10} ({joined});")
    regions_lines.append("}")

    region_solvers_lines = ["regionSolvers", "{"]
    for name, solver in solver_entries:
        region_solvers_lines.append(f"    {name:<12} {solver};")
    region_solvers_lines.append("}")

    return ("\n".join(regions_lines) + "\n", "\n".join(region_solvers_lines) + "\n")


def migrate_control_dict(text: str, case_dir: Path) -> str:
    blocks = regions_and_region_solvers_from_region_properties(case_dir / "constant/regionProperties")
    if not blocks:
        return text

    regions_block, region_solvers_block = blocks
    cleaned = re.sub(r"\n\s*regions\s*\{.*?\}\s*\n", "\n\n", text, flags=re.S)
    cleaned = re.sub(r"\n\s*regionSolvers\s*\{.*?\}\s*\n", "\n\n", cleaned, flags=re.S)
    insertion = "\n" + regions_block + "\n\n" + region_solvers_block + "\n"

    match = re.search(r"application\s+[^;]+;\s*\n", cleaned)
    if match:
        updated = cleaned[: match.end()] + insertion + cleaned[match.end() :]
    else:
        updated = regions_block + "\n" + region_solvers_block + "\n" + cleaned

    return updated


def ensure_fv_constraints(case_dir: Path, report: Report, apply: bool) -> None:
    air_system = case_dir / "system/air"
    fv_constraints = air_system / "fvConstraints"
    fv_options = air_system / "fvOptions"

    if not air_system.exists():
        return

    if not fv_constraints.exists():
        text = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system/air";
    object      fvConstraints;
}

limitp
{
    type        limitPressure;
    p           p_rgh;

    minFactor   0.8;
    maxFactor   1.2;
}
"""
        if apply:
            fv_constraints.write_text(text, encoding="utf-8")
        report.changed.append(str(fv_constraints))

    if fv_options.exists():
        opt_text = read_text(fv_options)
        if "limitTemperature" in opt_text:
            cons_text = read_text(fv_constraints) if fv_constraints.exists() else ""
            if "limitTemperature" not in cons_text:
                limit_match = re.search(r"\n\s*(\w+)\s*\{\s*type\s+limitTemperature;.*?\n\s*\}\s*", opt_text, flags=re.S)
                if limit_match:
                    limit_block = limit_match.group(0)
                    limit_block = re.sub(r"\bselectionMode\b", "select", limit_block)
                    cons_text = cons_text.rstrip() + "\n\n" + limit_block.strip() + "\n"
                    write_text(fv_constraints, cons_text, apply)
                    report.changed.append(str(fv_constraints))
                else:
                    report.review.append(f"manual migration needed for limitTemperature in {fv_options}")


def ensure_empty_fv_models(case_dir: Path, report: Report, apply: bool) -> None:
    air_constant = case_dir / "constant/air"
    if not air_constant.exists():
        return

    fv_models = air_constant / "fvModels"
    if fv_models.exists():
        return

    text = """FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant/air";
    object      fvModels;
}
"""
    if apply:
        fv_models.write_text(text, encoding="utf-8")
    report.changed.append(str(fv_models))


def maybe_rename_fvoptions(case_dir: Path, report: Report, apply: bool, keep: bool) -> None:
    if keep:
        return

    fv_options = case_dir / "system/air/fvOptions"
    backup = case_dir / "system/air/fvOptions.bak_of8"
    if not fv_options.exists():
        return

    if apply:
        if backup.exists():
            backup.unlink()
        shutil.move(str(fv_options), str(backup))
    report.changed.append(f"{fv_options} -> {backup}")


def system_view_factors_path(case_dir: Path, legacy_path: Path) -> Path | None:
    rel = legacy_path.relative_to(case_dir / "constant")
    parts = rel.parts

    if parts == ("viewFactorsDict",):
        return case_dir / "system/viewFactorsDict"

    if len(parts) == 2 and parts[1] == "viewFactorsDict":
        return case_dir / "system" / parts[0] / "viewFactorsDict"

    return None


def all_view_factors_targets(case_dir: Path, source_path: Path) -> list[Path]:
    targets: list[Path] = []

    if source_path.is_relative_to(case_dir / "constant"):
        rel = source_path.relative_to(case_dir / "constant")
        parts = rel.parts

        # OF12 faceAgglomerate in this build reads from system/viewFactorsDict,
        # while other radiation tools still read constant/<region>/viewFactorsDict.
        if parts == ("viewFactorsDict",):
            targets.append(case_dir / "system/viewFactorsDict")
        elif len(parts) == 2 and parts[1] == "viewFactorsDict":
            targets.append(case_dir / "system/viewFactorsDict")
            targets.append(case_dir / "system" / parts[0] / "viewFactorsDict")

    return targets


def relocate_view_factors_dicts(case_dir: Path, report: Report, apply: bool) -> None:
    constant_dir = case_dir / "constant"
    if not constant_dir.exists():
        return

    for source_path in sorted(constant_dir.rglob("viewFactorsDict")):
        if not source_path.is_file() or is_backup(source_path):
            continue

        targets = all_view_factors_targets(case_dir, source_path)
        if not targets:
            report.review.append(f"manual migration needed for viewFactorsDict at {source_path}")
            continue

        source_text = read_text(source_path)

        for target_path in targets:
            updated_text = source_text
            location = target_path.parent.relative_to(case_dir).as_posix()
            if re.search(r"(^\s*)location\s+\"[^\"]*\"\s*;", updated_text, flags=re.M):
                updated_text = re.sub(
                    r"(^\s*)location\s+\"[^\"]*\"\s*;",
                    rf'\1location    "{location}";',
                    updated_text,
                    flags=re.M,
                )

            if target_path.exists():
                existing_text = read_text(target_path)
                if existing_text != updated_text:
                    report.review.append(
                        f"viewFactorsDict conflict: {source_path} vs {target_path}"
                    )
                continue

            if apply:
                target_path.parent.mkdir(parents=True, exist_ok=True)
                target_path.write_text(updated_text, encoding="utf-8")

            report.changed.append(f"{source_path} -> {target_path}")


def relocate_model_property_dicts(case_dir: Path, report: Report, apply: bool) -> None:
    constant_dir = case_dir / "constant"
    air_dir = constant_dir / "air"
    if not constant_dir.exists() or not air_dir.exists():
        return

    if not (air_dir / "vegetationProperties").exists():
        report.review.append(f"missing vegetationProperties: {air_dir / 'vegetationProperties'}")
    if not (air_dir / "grassProperties").exists():
        report.review.append(f"missing grassProperties: {air_dir / 'grassProperties'}")


def translate_setset_line(line: str) -> dict[str, str] | None:
    stripped = line.strip()
    if not stripped or stripped.startswith("//") or stripped.startswith("#"):
        return None

    cell_set = re.fullmatch(
        r"cellSet\s+(\S+)\s+(new|add|delete)\s+boxToCell\s+(\([^)]*\))(\([^)]*\))",
        stripped,
    )
    if cell_set:
        name, action, box_min, box_max = cell_set.groups()
        return {
            "name": name,
            "type": "cellSet",
            "action": action,
            "source": "boxToCell",
            "box": f"{box_min} {box_max}",
        }

    cell_zone = re.fullmatch(
        r"cellZoneSet\s+(\S+)\s+(new|add|delete)\s+setToCellZone\s+(\S+)",
        stripped,
    )
    if cell_zone:
        name, action, set_name = cell_zone.groups()
        return {
            "name": name,
            "type": "cellZoneSet",
            "action": action,
            "source": "setToCellZone",
            "set": set_name,
        }

    return None


def topo_set_dict_text(location: str, actions: list[dict[str, str]]) -> str:
    lines = [
        "FoamFile",
        "{",
        "    format      ascii;",
        "    class       dictionary;",
        f'    location    "{location}";',
        "    object      topoSetDict;",
        "}",
        "",
        "actions",
        "(",
    ]

    for action in actions:
        lines.extend(
            [
                "    {",
                f'        name    {action["name"]};',
                f'        type    {action["type"]};',
                f'        action  {action["action"]};',
                f'        source  {action["source"]};',
            ]
        )
        if "box" in action:
            lines.append(f'        box     {action["box"]};')
        if "set" in action:
            lines.append(f'        set     {action["set"]};')
        lines.extend(["    }", ""])

    if lines[-1] == "":
        lines.pop()
    lines.extend([");", ""])
    return "\n".join(lines)


def migrate_setset_batches(case_dir: Path, report: Report, apply: bool) -> None:
    for batch_path in sorted((case_dir / "system").rglob("setset.batch")):
        if not batch_path.is_file() or is_backup(batch_path):
            continue

        region_dir = batch_path.parent
        target_path = region_dir / "topoSetDict"
        actions: list[dict[str, str]] = []
        unsupported: list[str] = []

        for raw_line in read_text(batch_path).splitlines():
            translated = translate_setset_line(raw_line)
            if translated is not None:
                actions.append(translated)
                continue

            if raw_line.strip():
                unsupported.append(raw_line.strip())

        if unsupported:
            report.review.append(
                f"manual migration needed for {batch_path}: unsupported setSet lines"
            )
            continue

        if not actions:
            report.review.append(f"manual migration needed for empty {batch_path}")
            continue

        location = target_path.parent.relative_to(case_dir).as_posix()
        updated_text = topo_set_dict_text(location, actions)

        if target_path.exists():
            existing_text = read_text(target_path)
            if existing_text != updated_text:
                report.review.append(f"topoSetDict conflict: {batch_path} vs {target_path}")
            continue

        if apply:
            target_path.write_text(updated_text, encoding="utf-8")
            backup_path = batch_path.with_suffix(batch_path.suffix + ".bak_of8")
            if backup_path.exists():
                backup_path.unlink()
            shutil.move(str(batch_path), str(backup_path))

        report.changed.append(f"{batch_path} -> {target_path}")


def review_scan(case_dir: Path, report: Report) -> None:
    patterns = {
        "legacy mapped patch keys": r"sampleMode|sampleRegion|samplePatch|type\s+mappedPatch|nearestPatchFace|nearestCell|nearestOnlyCell",
        "legacy pressure factors": r"pMaxFactor|pMinFactor",
        "legacy fvOptions header": r"object\s+fvOptions",
        "old RAS model key": r"RASModel|RASProperties",
        "old scripts": r"foamCleanPolyMesh|setSet",
        "old blockMesh substitution": r"\$:",
    }

    for path in all_files(case_dir):
        if path.name in RAW_TABLE_FILES or "processor" in path.parts:
            continue
        try:
            text = read_text(path)
        except UnicodeDecodeError:
            continue
        for label, pattern in patterns.items():
            if re.search(pattern, text):
                report.review.append(f"{label}: {path}")


def migrate_case(case_dir: Path, apply: bool, keep_fvoptions: bool) -> Report:
    report = Report()

    for path in all_files(case_dir):
        if path.name in RAW_TABLE_FILES or "processor" in path.parts:
            continue

        if path.name == "momentumTransport":
            update_file(path, migrate_momentum_transport, report, apply)
        elif path.name == "fvSolution":
            update_file(path, migrate_fvsolution, report, apply)
        elif path.name == "setFieldsDict":
            update_file(path, migrate_setfields, report, apply)
        elif path.name.startswith("blockMeshDict"):
            update_file(path, migrate_blockmesh, report, apply)
        elif path.name == "changeDictionaryDict":
            update_file(path, migrate_change_dictionary, report, apply)
        elif path.name == "controlDict" and path.parent.name == "system":
            update_file(path, lambda text: migrate_control_dict(text, case_dir), report, apply)

    ensure_fv_constraints(case_dir, report, apply)
    ensure_empty_fv_models(case_dir, report, apply)
    maybe_rename_fvoptions(case_dir, report, apply, keep_fvoptions)
    relocate_view_factors_dicts(case_dir, report, apply)
    relocate_model_property_dicts(case_dir, report, apply)
    migrate_setset_batches(case_dir, report, apply)
    migrate_region_initial_fields(case_dir, report, apply)
    ensure_run_scripts(case_dir, report, apply)
    review_scan(case_dir, report)

    return report


def main() -> int:
    args = parse_args()
    case_dir = args.case.resolve()
    if not case_dir.is_dir():
        raise SystemExit(f"case directory does not exist: {case_dir}")

    report = migrate_case(case_dir, args.apply, args.keep_fvoptions)
    mode = "applied" if args.apply else "dry-run"
    print(f"OF12 migration {mode}: {case_dir}")

    if report.changed:
        print("\nChanges:")
        for item in report.changed:
            prefix = "changed" if args.apply else "would change"
            print(f"  {prefix}: {item}")
    else:
        print("\nChanges: none")

    if report.review:
        print("\nManual review:")
        for item in sorted(set(report.review)):
            print(f"  {item}")
    else:
        print("\nManual review: no known legacy patterns remain")

    if not args.apply:
        print("\nDry-run only. Re-run with --apply to write changes.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

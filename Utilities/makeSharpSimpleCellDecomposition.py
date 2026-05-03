#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Erzeugt constant/cellDecomposition aus constant/polyMesh
für eine scharfe, blockartige Zerlegung wie simpleCoeffs n (Nx Ny Nz),
aber als manual decomposition.

Ziel:
- keine cellZones nötig
- beliebig z.B. (2 1 2), (10 1 20), ...
- Schnitte werden global in x/y/z gesetzt
- Schnitte werden auf komplette Zellschichten geschnappt
- Priorität: gerade/blockartige Schnitte, nicht exakt gleiche Zellanzahl

Aufruf:
    python3 makeSharpSimpleDecomposition.py --n 2 1 2
    python3 makeSharpSimpleDecomposition.py --n 10 1 20 /pfad/zum/case

Danach in system/decomposeParDict:
    numberOfSubdomains 4;
    method manual;

    manualCoeffs
    {
        dataFile "cellDecomposition";
    }
"""

from __future__ import annotations

import argparse
import bisect
import re
import sys
from collections import Counter
from pathlib import Path
from typing import List, Tuple, Set


OUTPUT_NAME = "cellDecomposition"


# -----------------------------------------------------------------------------
# Allgemein
# -----------------------------------------------------------------------------

def die(msg: str, code: int = 1) -> None:
    print(f"FEHLER: {msg}", file=sys.stderr)
    sys.exit(code)


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//.*", "", text)
    return text


def find_matching(text: str, start: int, open_char: str, close_char: str) -> int:
    if start < 0 or start >= len(text) or text[start] != open_char:
        raise ValueError(f"Position {start} zeigt nicht auf '{open_char}'.")

    depth = 0
    for i in range(start, len(text)):
        if text[i] == open_char:
            depth += 1
        elif text[i] == close_char:
            depth -= 1
            if depth == 0:
                return i

    raise ValueError(f"Kein passendes '{close_char}' gefunden.")


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        die(f"Datei nicht gefunden: {path}")
    except Exception as e:
        die(f"Konnte Datei nicht lesen: {path}\n{e}")


def remove_foamfile_header(text: str) -> str:
    text = strip_comments(text)

    if re.search(r"\bformat\s+binary\s*;", text):
        die(
            "Dieses Skript erwartet ascii-polyMesh-Dateien. "
            "Konvertiere ggf. zuerst mit foamFormatConvert oder schreibe das Mesh als ascii."
        )

    m = re.search(r"\bFoamFile\b", text)
    if not m:
        return text

    brace_open = text.find("{", m.end())
    if brace_open == -1:
        return text

    brace_close = find_matching(text, brace_open, "{", "}")
    return text[brace_close + 1:]


def read_foam_list_blob(path: Path) -> Tuple[int, str]:
    text = remove_foamfile_header(read_text(path))

    m = re.search(r"\b(\d+)\s*\(", text, flags=re.S)
    if not m:
        die(f"Konnte OpenFOAM-Liste nicht lesen: {path}")

    n_declared = int(m.group(1))
    open_pos = text.find("(", m.end() - 1)
    close_pos = find_matching(text, open_pos, "(", ")")

    return n_declared, text[open_pos + 1:close_pos]


# -----------------------------------------------------------------------------
# polyMesh lesen
# -----------------------------------------------------------------------------

def read_points(path: Path) -> List[Tuple[float, float, float]]:
    n_declared, blob = read_foam_list_blob(path)

    vec_re = re.compile(
        r"\(\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s+"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s+"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
        r"\s*\)"
    )

    points = [
        (float(a), float(b), float(c))
        for a, b, c in vec_re.findall(blob)
    ]

    if len(points) != n_declared:
        die(f"points: deklariert {n_declared}, gelesen {len(points)}.")

    return points


def read_label_list(path: Path) -> List[int]:
    n_declared, blob = read_foam_list_blob(path)
    values = [int(x) for x in re.findall(r"-?\d+", blob)]

    if len(values) != n_declared:
        die(f"{path.name}: deklariert {n_declared}, gelesen {len(values)}.")

    return values


def read_faces(path: Path) -> List[List[int]]:
    n_declared, blob = read_foam_list_blob(path)

    faces: List[List[int]] = []
    face_re = re.compile(r"(\d+)\s*\(([^()]*)\)", flags=re.S)

    for n_str, ids_blob in face_re.findall(blob):
        n_face_points = int(n_str)
        ids = [int(x) for x in re.findall(r"-?\d+", ids_blob)]

        if len(ids) != n_face_points:
            die(
                f"faces: Face erwartet {n_face_points} Punkte, "
                f"gelesen {len(ids)}."
            )

        faces.append(ids)

    if len(faces) != n_declared:
        die(f"faces: deklariert {n_declared}, gelesen {len(faces)}.")

    return faces


def build_cell_points(
    faces: List[List[int]],
    owner: List[int],
    neighbour: List[int],
) -> List[Set[int]]:
    if len(owner) != len(faces):
        die(f"owner hat {len(owner)} Einträge, faces hat {len(faces)} Einträge.")

    max_cell = max(owner + neighbour) if neighbour else max(owner)
    n_cells = max_cell + 1

    cell_points: List[Set[int]] = [set() for _ in range(n_cells)]

    for face_i, pts in enumerate(faces):
        own = owner[face_i]
        cell_points[own].update(pts)

        # In OpenFOAM gehören die neighbour-Einträge zu den ersten nInternalFaces.
        if face_i < len(neighbour):
            nei = neighbour[face_i]
            cell_points[nei].update(pts)

    for cell_i, pts in enumerate(cell_points):
        if not pts:
            die(f"Zelle {cell_i} hat keine Punkte bekommen.")

    return cell_points


def compute_cell_centres(
    points: List[Tuple[float, float, float]],
    cell_points: List[Set[int]],
) -> List[Tuple[float, float, float]]:
    centres: List[Tuple[float, float, float]] = []

    for cell_i, ids in enumerate(cell_points):
        try:
            xs = [points[p][0] for p in ids]
            ys = [points[p][1] for p in ids]
            zs = [points[p][2] for p in ids]
        except IndexError:
            die(f"Zelle {cell_i} referenziert einen ungültigen Punkt.")

        n = float(len(ids))
        centres.append((sum(xs) / n, sum(ys) / n, sum(zs) / n))

    return centres


# -----------------------------------------------------------------------------
# Achsen-Schnitte
# -----------------------------------------------------------------------------

def cluster_axis_values(values: List[float], tol_abs: float) -> Tuple[List[float], List[int]]:
    """
    Gruppiert fast gleiche Koordinaten zu Zellschichten.

    Rückgabe:
    - levels: sortierte repräsentative Koordinate je Schicht
    - ranks[cell] = Schichtindex der Zelle
    """
    pairs = sorted((v, i) for i, v in enumerate(values))

    levels: List[float] = []
    ranks = [-1] * len(values)

    current_indices: List[int] = []
    current_values: List[float] = []

    def flush() -> None:
        if not current_indices:
            return
        level_index = len(levels)
        level_value = sum(current_values) / len(current_values)
        levels.append(level_value)
        for _, cell_i in current_indices:
            ranks[cell_i] = level_index

    for v, cell_i in pairs:
        if not current_values:
            current_values = [v]
            current_indices = [(v, cell_i)]
            continue

        current_mean = sum(current_values) / len(current_values)

        if abs(v - current_mean) <= tol_abs:
            current_values.append(v)
            current_indices.append((v, cell_i))
        else:
            flush()
            current_values = [v]
            current_indices = [(v, cell_i)]

    flush()

    if any(r < 0 for r in ranks):
        die("Interner Fehler beim Clustering der Achsenwerte.")

    return levels, ranks


def make_rank_cuts(
    levels: List[float],
    n_parts: int,
    domain_min: float,
    domain_max: float,
    strategy: str,
) -> List[int]:
    """
    Erzeugt Schnittindizes in der sortierten Level-Liste.

    cuts = [0, ..., nLevels]
    Ein Level rank gehört zu part p, wenn cuts[p] <= rank < cuts[p+1].
    """
    n_levels = len(levels)

    if n_parts < 1:
        die("n_parts muss >= 1 sein.")

    if n_parts > n_levels:
        die(
            f"Es sollen {n_parts} Teile erzeugt werden, aber es gibt nur "
            f"{n_levels} Zellschichten in dieser Richtung."
        )

    if n_parts == 1:
        return [0, n_levels]

    cuts = [0]
    last = 0

    for p in range(1, n_parts):
        if strategy == "equal-levels":
            idx = round(p * n_levels / n_parts)
        elif strategy == "geometric":
            ideal = domain_min + (domain_max - domain_min) * p / n_parts
            idx = bisect.bisect_left(levels, ideal)
        else:
            die(f"Unbekannte strategy: {strategy}")

        # Nichtleere Teile erzwingen.
        min_idx = last + 1
        max_idx = n_levels - (n_parts - p)
        idx = max(min_idx, min(idx, max_idx))

        cuts.append(idx)
        last = idx

    cuts.append(n_levels)
    return cuts


def rank_to_part(rank: int, cuts: List[int]) -> int:
    return bisect.bisect_right(cuts, rank) - 1


def build_axis_parts(
    centres: List[Tuple[float, float, float]],
    points: List[Tuple[float, float, float]],
    axis: int,
    n_parts: int,
    tol_rel: float,
    tol_abs_min: float,
    strategy: str,
) -> Tuple[List[int], List[float], List[int]]:
    centre_values = [c[axis] for c in centres]
    point_values = [p[axis] for p in points]

    domain_min = min(point_values)
    domain_max = max(point_values)
    span = max(domain_max - domain_min, 1.0)

    tol_abs = max(tol_abs_min, tol_rel * span)

    levels, ranks = cluster_axis_values(centre_values, tol_abs)
    cuts = make_rank_cuts(levels, n_parts, domain_min, domain_max, strategy)

    parts = [rank_to_part(r, cuts) for r in ranks]

    return parts, levels, cuts


# -----------------------------------------------------------------------------
# Schreiben
# -----------------------------------------------------------------------------

def write_cell_decomposition(path: Path, decomp: List[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        f.write(
            "FoamFile\n"
            "{\n"
            "    version     2.0;\n"
            "    format      ascii;\n"
            "    class       labelList;\n"
            f"    object      {OUTPUT_NAME};\n"
            "}\n"
            "\n"
        )

        f.write(f"{len(decomp)}\n")
        f.write("(\n")
        for p in decomp:
            f.write(f"{p}\n")
        f.write(")\n")


def approx_cut_positions(levels: List[float], cuts: List[int]) -> List[float]:
    result: List[float] = []

    for idx in cuts[1:-1]:
        if idx <= 0:
            result.append(levels[0])
        elif idx >= len(levels):
            result.append(levels[-1])
        else:
            result.append(0.5 * (levels[idx - 1] + levels[idx]))

    return result


def print_summary(
    case: Path,
    n_xyz: Tuple[int, int, int],
    decomp: List[int],
    axis_info: List[Tuple[str, List[float], List[int]]],
) -> None:
    nx, ny, nz = n_xyz
    n_subdomains = nx * ny * nz
    counts = Counter(decomp)

    print(f"geschrieben        : {case / 'constant' / OUTPUT_NAME}")
    print(f"nCells             : {len(decomp)}")
    print(f"n                  : ({nx} {ny} {nz})")
    print(f"numberOfSubdomains : {n_subdomains}")
    print()

    print("Achsen-Schnitte:")
    for name, levels, cuts in axis_info:
        cut_pos = approx_cut_positions(levels, cuts)
        print(f"  {name}: {len(levels)} Zellschichten")
        if cut_pos:
            print("     Schnittpositionen ungefähr:", " ".join(f"{x:.12g}" for x in cut_pos))
        else:
            print("     keine Schnitte")

    print()
    print("Processor -> Anzahl Zellen:")
    for p in range(n_subdomains):
        print(f"  {p:>5} -> {counts.get(p, 0)}")

    empty = [p for p in range(n_subdomains) if counts.get(p, 0) == 0]
    if empty:
        print()
        print("WARNUNG: Leere Prozessoren:", empty)

    print()
    print("system/decomposeParDict:")
    print()
    print(f"numberOfSubdomains {n_subdomains};")
    print("method manual;")
    print()
    print("manualCoeffs")
    print("{")
    print(f'    dataFile "{OUTPUT_NAME}";')
    print("}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Erzeugt eine scharfe simple-artige manual cellDecomposition."
    )

    parser.add_argument(
        "case",
        nargs="?",
        default=".",
        help="Pfad zum OpenFOAM-Case, Standard: aktuelles Verzeichnis",
    )

    parser.add_argument(
        "--n",
        nargs=3,
        type=int,
        metavar=("NX", "NY", "NZ"),
        required=True,
        help="Zerlegung wie simpleCoeffs n, z.B. --n 2 1 2",
    )

    parser.add_argument(
        "--strategy",
        choices=["geometric", "equal-levels"],
        default="geometric",
        help=(
            "geometric: ideale Schnitte wie simple nach Domain-Länge, "
            "aber auf Zellschichten geschnappt. "
            "equal-levels: gleiche Anzahl Zellschichten je Richtung."
        ),
    )

    parser.add_argument(
        "--axisTolRel",
        type=float,
        default=1e-8,
        help="Relative Toleranz zum Gruppieren fast gleicher Zellschicht-Koordinaten.",
    )

    parser.add_argument(
        "--axisTolAbs",
        type=float,
        default=1e-14,
        help="Minimale absolute Toleranz zum Gruppieren fast gleicher Koordinaten.",
    )

    args = parser.parse_args()

    nx, ny, nz = args.n
    if nx < 1 or ny < 1 or nz < 1:
        die("--n muss positive Werte haben.")

    case = Path(args.case).resolve()
    poly = case / "constant" / "polyMesh"

    points_file = poly / "points"
    faces_file = poly / "faces"
    owner_file = poly / "owner"
    neighbour_file = poly / "neighbour"

    if not poly.exists():
        die(f"polyMesh-Verzeichnis nicht gefunden: {poly}")

    points = read_points(points_file)
    faces = read_faces(faces_file)
    owner = read_label_list(owner_file)

    if neighbour_file.exists():
        neighbour = read_label_list(neighbour_file)
    else:
        neighbour = []

    cell_points = build_cell_points(faces, owner, neighbour)
    centres = compute_cell_centres(points, cell_points)

    x_parts, x_levels, x_cuts = build_axis_parts(
        centres, points, 0, nx, args.axisTolRel, args.axisTolAbs, args.strategy
    )
    y_parts, y_levels, y_cuts = build_axis_parts(
        centres, points, 1, ny, args.axisTolRel, args.axisTolAbs, args.strategy
    )
    z_parts, z_levels, z_cuts = build_axis_parts(
        centres, points, 2, nz, args.axisTolRel, args.axisTolAbs, args.strategy
    )

    decomp: List[int] = []
    for cell_i in range(len(centres)):
        ix = x_parts[cell_i]
        iy = y_parts[cell_i]
        iz = z_parts[cell_i]

        # OpenFOAM-simple-artige Nummerierung:
        # x läuft am schnellsten, dann y, dann z.
        proc = ix + nx * (iy + ny * iz)
        decomp.append(proc)

    expected = nx * ny * nz
    if max(decomp) >= expected:
        die("Interner Fehler: Prozessor-ID außerhalb des erwarteten Bereichs.")

    out_file = case / "constant" / OUTPUT_NAME
    write_cell_decomposition(out_file, decomp)

    print_summary(
        case,
        (nx, ny, nz),
        decomp,
        [
            ("x", x_levels, x_cuts),
            ("y", y_levels, y_cuts),
            ("z", z_levels, z_cuts),
        ],
    )


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Erzeugt constant/cellDecomposition aus constant/polyMesh/cellZones
für method manual in decomposePar.

Standardverhalten:
- liest alle cellZones aus constant/polyMesh/cellZones
- liest die Gesamtzahl der Zellen aus constant/polyMesh/cells
- weist allen Zellen einer Zone dieselbe Prozessor-ID zu
- schreibt constant/cellDecomposition

Wichtige Annahme:
- Jede Zelle, die verteilt werden soll, muss genau einer Zone zugeordnet sein.
- Wenn Zellen keiner Zone zugeordnet sind, bricht das Script standardmäßig ab.
  Das kann unten über ASSIGN_UNZONED_TO geändert werden.

Aufruf:
    python3 makeCellDecomposition.py

Optional:
    python3 makeCellDecomposition.py /pfad/zum/case
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional


# -----------------------------------------------------------------------------
# KONFIGURATION
# -----------------------------------------------------------------------------

# Falls None:
#   Jede Zone bekommt automatisch eine eigene Prozessor-ID in der Reihenfolge
#   des cellZones-Files: erste Zone -> 0, zweite Zone -> 1, ...
#
# Beispiel für manuelles Mapping:
# write the zone name to mapped processor
ZONE_TO_PROC = {
    "zR4S1": 0,
    "zR4S2": 1,
    "zR4S3": 2,
    "zR4S4": 3,
    "zR3S1": 4,
    "zR3S2": 5,
    "zR3S3": 6,
    "zR3S4": 7,
    "zR2S1": 4,
    "zR2S2": 5,
    "zR2S3": 6,
    "zR2S4": 7,
    "zR1S1": 8,
}
# for automatic mapping
#ZONE_TO_PROC: Optional[Dict[str, int]] = None   

# Wenn None:
#   Script bricht ab, sobald es Zellen findet, die in keiner cellZone liegen.
#
# Wenn z.B. 0:
#   alle unzugeordneten Zellen werden auf Prozessor 0 gelegt.
ASSIGN_UNZONED_TO: Optional[int] = None

# Optionaler Plausibilitätscheck:
# Wenn gesetzt, muss max(proc)+1 genau diesem Wert entsprechen.
# Sonst bricht das Script ab.
EXPECTED_NUMBER_OF_SUBDOMAINS: Optional[int] = 9 # oder None für keinen Check

# Name der Ausgabedatei
OUTPUT_NAME = "cellDecomposition"


# -----------------------------------------------------------------------------
# HILFSFUNKTIONEN
# -----------------------------------------------------------------------------

def die(msg: str, code: int = 1) -> None:
    print(f"FEHLER: {msg}", file=sys.stderr)
    sys.exit(code)


def strip_comments(text: str) -> str:
    # Block-Kommentare
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    # Zeilen-Kommentare
    text = re.sub(r"//.*", "", text)
    return text


def find_matching(text: str, start: int, open_char: str, close_char: str) -> int:
    if start < 0 or start >= len(text) or text[start] != open_char:
        raise ValueError(
            f"find_matching: Position {start} zeigt nicht auf '{open_char}'."
        )

    depth = 0
    for i in range(start, len(text)):
        c = text[i]
        if c == open_char:
            depth += 1
        elif c == close_char:
            depth -= 1
            if depth == 0:
                return i

    raise ValueError(
        f"Kein passendes schließendes '{close_char}' für '{open_char}' bei {start} gefunden."
    )


def remove_foamfile_header(text: str) -> str:
    text = strip_comments(text)

    m = re.search(r"\bFoamFile\b", text)
    if not m:
        return text

    brace_open = text.find("{", m.end())
    if brace_open == -1:
        return text

    brace_close = find_matching(text, brace_open, "{", "}")
    return text[brace_close + 1:]


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        die(f"Datei nicht gefunden: {path}")
    except Exception as e:
        die(f"Konnte Datei nicht lesen: {path}\n{e}")


# -----------------------------------------------------------------------------
# POLYMESH EINLESEN
# -----------------------------------------------------------------------------

def read_list_length(path: Path) -> int:
    """
    Liest die Länge einer OpenFOAM-Liste:
        N
        (
            ...
        )
    z.B. aus cells oder cellLevel.
    """
    text = remove_foamfile_header(read_text(path))
    m = re.search(r"\b(\d+)\s*\(", text, flags=re.S)
    if not m:
        die(f"Konnte Listenlänge in {path} nicht lesen.")
    n = int(m.group(1))
    if n <= 0:
        die(f"Ungültige Listenlänge in {path}: {n}")
    return n


def read_ncells(poly_mesh: Path, zones: List[Tuple[str, List[int]]]) -> int:
    """
    Bestimmt die Zellanzahl robust:
    1) cells      (falls vorhanden)
    2) cellLevel  (falls vorhanden)
    3) max(cellZone-ID)+1 als Fallback
    """
    cells_file = poly_mesh / "cells"
    if cells_file.exists():
        return read_list_length(cells_file)

    cell_level_file = poly_mesh / "cellLevel"
    if cell_level_file.exists():
        return read_list_length(cell_level_file)

    # Fallback: nur aus den cellZones
    max_id = -1
    for zone_name, cell_ids in zones:
        if cell_ids:
            max_id = max(max_id, max(cell_ids))

    if max_id < 0:
        die(
            "Konnte nCells weder aus cells noch aus cellLevel bestimmen, "
            "und die cellZones sind leer."
        )

    return max_id + 1


def parse_cellzones(cellzones_file: Path) -> List[Tuple[str, List[int]]]:
    """
    Liest constant/polyMesh/cellZones.
    Erwartet pro Zone mindestens einen Block:
        zoneName
        {
            ...
            cellLabels List<label> N ( ... );
            ...
        }
    """
    text = remove_foamfile_header(read_text(cellzones_file))

    # Äußerer Listenblock: z.B.
    # 13
    # (
    #   zoneA { ... }
    #   zoneB { ... }
    # )
    m = re.search(r"\b(\d+)\s*\(", text, flags=re.S)
    if not m:
        die(f"Konnte äußeren Zonenblock in {cellzones_file} nicht finden.")

    n_zones_declared = int(m.group(1))
    list_open = text.find("(", m.end() - 1)
    if list_open == -1:
        die(f"Konnte öffnende Klammer des Zonenblocks in {cellzones_file} nicht finden.")

    list_close = find_matching(text, list_open, "(", ")")
    zones_blob = text[list_open + 1:list_close]

    zones: List[Tuple[str, List[int]]] = []

    i = 0
    n = len(zones_blob)

    while i < n:
        # Leerraum überspringen
        while i < n and zones_blob[i].isspace():
            i += 1
        if i >= n:
            break

        # Zonenname
        name_match = re.match(r"([A-Za-z_]\w*)", zones_blob[i:])
        if not name_match:
            snippet = zones_blob[i:i + 80].replace("\n", " ")
            die(
                f"Konnte Zonenname in cellZones nicht lesen bei:\n{snippet}"
            )

        zone_name = name_match.group(1)
        i += name_match.end()

        # Leerraum bis '{'
        while i < n and zones_blob[i].isspace():
            i += 1

        if i >= n or zones_blob[i] != "{":
            die(f"Nach Zonenname '{zone_name}' wurde '{{' erwartet.")

        body_open = i
        body_close = find_matching(zones_blob, body_open, "{", "}")
        zone_body = zones_blob[body_open + 1:body_close]
        i = body_close + 1

        # cellLabels finden
        cm = re.search(
            r"cellLabels\s+List<label>\s+(\d+)\s*\(",
            zone_body,
            flags=re.S
        )
        if not cm:
            die(
                f"In Zone '{zone_name}' wurde kein 'cellLabels List<label> ... (...)' gefunden."
            )

        n_labels_declared = int(cm.group(1))
        labels_open = zone_body.find("(", cm.end() - 1)
        if labels_open == -1:
            die(f"In Zone '{zone_name}' konnte die cellLabels-Liste nicht geöffnet werden.")

        labels_close = find_matching(zone_body, labels_open, "(", ")")
        labels_blob = zone_body[labels_open + 1:labels_close]

        labels = [int(x) for x in re.findall(r"-?\d+", labels_blob)]

        if len(labels) != n_labels_declared:
            die(
                f"Zone '{zone_name}': erwartet {n_labels_declared} Zellen, "
                f"gefunden {len(labels)}."
            )

        zones.append((zone_name, labels))

    if len(zones) != n_zones_declared:
        die(
            f"cellZones: deklariert sind {n_zones_declared} Zonen, "
            f"gelesen wurden {len(zones)}."
        )

    return zones


# -----------------------------------------------------------------------------
# DECOMPOSITION ERZEUGEN
# -----------------------------------------------------------------------------

def build_zone_to_proc(zones: List[Tuple[str, List[int]]]) -> Dict[str, int]:
    zone_names = [name for name, _ in zones]

    if ZONE_TO_PROC is None:
        return {name: i for i, name in enumerate(zone_names)}

    # Prüfen, ob alle Zonen ein Mapping haben
    missing = [name for name in zone_names if name not in ZONE_TO_PROC]
    if missing:
        die(
            "Für folgende Zonen fehlt ein Prozessor-Mapping in ZONE_TO_PROC:\n"
            + "\n".join(f"  {z}" for z in missing)
        )

    # Prüfen, ob Mapping negative IDs enthält
    bad = {k: v for k, v in ZONE_TO_PROC.items() if not isinstance(v, int) or v < 0}
    if bad:
        lines = [f"  {k}: {v}" for k, v in bad.items()]
        die("Ungültige Prozessor-IDs in ZONE_TO_PROC:\n" + "\n".join(lines))

    return dict(ZONE_TO_PROC)


def make_decomposition(
    n_cells: int,
    zones: List[Tuple[str, List[int]]],
    zone_to_proc: Dict[str, int],
) -> List[int]:
    decomp = [-1] * n_cells

    for zone_name, cell_ids in zones:
        proc = zone_to_proc[zone_name]

        for cell_id in cell_ids:
            if cell_id < 0 or cell_id >= n_cells:
                die(
                    f"Zone '{zone_name}' enthält ungültige Zelle {cell_id} "
                    f"(gültig: 0..{n_cells - 1})."
                )

            if decomp[cell_id] != -1:
                die(
                    f"Zelle {cell_id} liegt in mehreren Zonen: "
                    f"mindestens in Prozessor {decomp[cell_id]} und zusätzlich in Zone '{zone_name}'."
                )

            decomp[cell_id] = proc

    unassigned = [i for i, p in enumerate(decomp) if p == -1]

    if unassigned:
        if ASSIGN_UNZONED_TO is None:
            preview = ", ".join(str(x) for x in unassigned[:20])
            die(
                f"{len(unassigned)} Zellen sind keiner cellZone zugeordnet.\n"
                f"Beispiele: {preview}\n"
                f"Setze ASSIGN_UNZONED_TO auf eine Prozessor-ID, falls diese Zellen "
                f"bewusst auf einen festen Prozessor gelegt werden sollen."
            )

        if not isinstance(ASSIGN_UNZONED_TO, int) or ASSIGN_UNZONED_TO < 0:
            die(f"ASSIGN_UNZONED_TO ist ungültig: {ASSIGN_UNZONED_TO}")

        for cell_id in unassigned:
            decomp[cell_id] = ASSIGN_UNZONED_TO

    return decomp


def write_cell_decomposition(out_file: Path, decomp: List[int]) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        with out_file.open("w", encoding="utf-8") as f:
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
    except Exception as e:
        die(f"Konnte {out_file} nicht schreiben.\n{e}")


# -----------------------------------------------------------------------------
# AUSGABE
# -----------------------------------------------------------------------------

def print_summary(
    n_cells: int,
    zones: List[Tuple[str, List[int]]],
    zone_to_proc: Dict[str, int],
    decomp: List[int],
    out_file: Path,
) -> None:
    proc_counts = Counter(decomp)
    n_subdomains = max(decomp) + 1 if decomp else 0

    if EXPECTED_NUMBER_OF_SUBDOMAINS is not None:
        if n_subdomains != EXPECTED_NUMBER_OF_SUBDOMAINS:
            die(
                f"Erzeugte numberOfSubdomains={n_subdomains}, erwartet war "
                f"{EXPECTED_NUMBER_OF_SUBDOMAINS}. "
                f"Prüfe ZONE_TO_PROC oder EXPECTED_NUMBER_OF_SUBDOMAINS."
            )

    print(f"geschrieben           : {out_file}")
    print(f"nCells                : {n_cells}")
    print(f"nZones                : {len(zones)}")
    print(f"numberOfSubdomains    : {n_subdomains}")
    print()

    print("Zone -> Prozessor:")
    max_name = max(len(name) for name, _ in zones) if zones else 10
    for zone_name, cell_ids in zones:
        print(
            f"  {zone_name:<{max_name}} -> {zone_to_proc[zone_name]:>3}   "
            f"({len(cell_ids)} cells)"
        )

    print()
    print("Processor -> Anzahl Zellen:")
    for proc in range(n_subdomains):
        print(f"  {proc:>3} -> {proc_counts.get(proc, 0)}")

    print()
    print("Trage in system/decomposeParDict ein:")
    print()
    print("numberOfSubdomains ", n_subdomains, ";", sep="")
    print("method manual;")
    print()
    print("manualCoeffs")
    print("{")
    print(f'    dataFile "{OUTPUT_NAME}";')
    print("}")


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Erzeugt constant/cellDecomposition aus constant/polyMesh/cellZones."
    )
    parser.add_argument(
        "case",
        nargs="?",
        default=".",
        help="Pfad zum OpenFOAM-Case (Standard: aktuelles Verzeichnis)",
    )
    args = parser.parse_args()

    case = Path(args.case).resolve()

    poly_mesh = case / "constant" / "polyMesh"
    cells_file = poly_mesh / "cells"
    cellzones_file = poly_mesh / "cellZones"
    out_file = case / "constant" / OUTPUT_NAME

    if not poly_mesh.exists():
        die(f"polyMesh-Verzeichnis nicht gefunden: {poly_mesh}")

    zones = parse_cellzones(cellzones_file)
    n_cells = read_ncells(poly_mesh, zones)
    zone_to_proc = build_zone_to_proc(zones)
    decomp = make_decomposition(n_cells, zones, zone_to_proc)
    write_cell_decomposition(out_file, decomp)
    print_summary(n_cells, zones, zone_to_proc, decomp, out_file)


if __name__ == "__main__":
    main()
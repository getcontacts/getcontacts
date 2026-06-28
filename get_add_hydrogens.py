"""
get_add_hydrogens.py — Add hydrogens to a protein/nucleic-acid structure file.

Pure-Python reimplementation of the Richardson-lab "reduce" utility.
Reads/writes PDB format. No external dependencies beyond the standard library.

ORIGINS
-------
This file is an independent Python reimplementation. The algorithms, geometry
tables, and scoring heuristics are derived from studying two upstream projects:

  1. reduce (C++ original)
     Richardson Lab / J. Michael Word
     https://github.com/phenix-project/reduce/tree/master/reduce_src

  2. mmtbx.reduce (partial Python reimplementation in CCTBX)
     The Regents of the University of California, through Lawrence Berkeley
     National Laboratory
     https://github.com/cctbx/cctbx_project/tree/master/mmtbx/reduce

No source code was copied from either project. All geometry, scoring, and
optimization logic was rewritten from scratch in Python based on the published
algorithms and documentation.

LICENSE NOTICES
---------------
Neither upstream project uses a copyleft (GPL-style) license, so there is no
obligation to apply their terms to this file. Their notices are reproduced here
as a courtesy and in the spirit of the original requests for attribution.

  reduce (phenix-project):
    Copyright (C) 1999 J. Michael Word
    Copyright (C) 2020 ReliaSolve
    Copyright (c) 1989, 1992 The Regents of the University of California.
    All rights reserved.

    NOTICE: This is free software and the source code is freely available.
    You are free to redistribute or modify under the conditions that
    (1) this notice is not removed or modified in any way and
    (2) any modified versions of the program are also available for free.
    ** Absolutely no Warranty **

    Redistribution and use in source and binary forms are permitted provided
    that the above copyright notice and this paragraph are duplicated in all
    such forms and that any documentation, advertising materials, and other
    materials related to such distribution and use acknowledge that the software
    was developed by the University of California, San Francisco. The name of
    the University may not be used to endorse or promote products derived from
    this software without specific prior written permission.

  CCTBX / mmtbx.reduce:
    cctbx Copyright (c) 2006-2026, The Regents of the University of California,
    through Lawrence Berkeley National Laboratory. All rights reserved.
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that: (1) redistributions of source
    code retain the above copyright notice, this list of conditions and the
    disclaimer; (2) redistributions in binary form reproduce the above copyright
    notice in the documentation and/or other materials provided; (3) neither the
    name of the University of California, Lawrence Berkeley National Laboratory,
    nor U.S. Dept. of Energy nor the names of contributors may be used to
    endorse or promote products derived from this software without specific prior
    written permission. THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF
    ANY KIND.

API
---
    from get_add_hydrogens import add_hydrogens
    add_hydrogens("input.pdb", "output.pdb")

CLI
---
    python get_add_hydrogens.py input.pdb output.pdb [--no-flip] [--no-optimize]
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Bond lengths (Å), X-ray values (engh-huber / standard small-molecule)
BL_CH_ALI  = 0.970   # aliphatic C-H
BL_CH_ARO  = 0.930   # aromatic C-H
BL_NH      = 0.860   # N-H (amine/amide)
BL_OH      = 0.840   # O-H
BL_SH      = 1.200   # S-H
BL_CH2_ALI = 0.970   # methylene C-H (same as aliphatic)

# Bond angles (degrees)
TETRA = 109.471       # perfect tetrahedral
PLANAR = 120.0        # sp2 planar
LINEAR = 180.0

# VDW radii (Å) — from Bondi / reduce parameter sets
VDW = {
    "H": 1.17, "C": 1.70, "N": 1.55, "O": 1.52,
    "S": 1.80, "P": 1.80, "F": 1.47, "CL": 1.75,
    "BR": 1.85, "I": 1.98, "SE": 1.90,
}

# Clash threshold: fraction of sum-of-VDW below which it's a clash
CLASH_OVERLAP_CUTOFF = 0.40   # Å — contacts worse than this count as clashes

# H-bond geometry
HBOND_DIST_CUTOFF = 3.50      # donor-acceptor distance
HBOND_ANGLE_CUTOFF = 90.0     # min donor-H…acceptor angle (degrees)

# ---------------------------------------------------------------------------
# PDB I/O
# ---------------------------------------------------------------------------

@dataclass
class Atom:
    record:   str   = "ATOM  "   # "ATOM  " or "HETATM"
    serial:   int   = 0
    name:     str   = ""         # e.g. " CA " (4-char PDB style)
    altloc:   str   = " "
    resname:  str   = ""
    chain:    str   = "A"
    resseq:   int   = 0
    icode:    str   = " "
    x:        float = 0.0
    y:        float = 0.0
    z:        float = 0.0
    occupancy:float = 1.0
    bfactor:  float = 0.0
    element:  str   = ""
    charge:   str   = "  "

    @property
    def pos(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def name_stripped(self) -> str:
        return self.name.strip()

    def to_pdb_line(self) -> str:
        elem = self.element.rjust(2) if self.element else "  "
        return (
            f"{self.record:<6}"
            f"{self.serial:5d} "
            f"{self.name:<4}"
            f"{self.altloc}"
            f"{self.resname:<3} "
            f"{self.chain}"
            f"{self.resseq:4d}"
            f"{self.icode}   "
            f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"
            f"{self.occupancy:6.2f}{self.bfactor:6.2f}          "
            f"{elem}"
            f"{self.charge}\n"
        )


def parse_pdb(lines: List[str]) -> List[Atom]:
    atoms: List[Atom] = []
    for line in lines:
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        try:
            a = Atom(
                record   = rec,
                serial   = int(line[6:11]),
                name     = line[12:16],
                altloc   = line[16],
                resname  = line[17:20].strip(),
                chain    = line[21],
                resseq   = int(line[22:26]),
                icode    = line[26],
                x        = float(line[30:38]),
                y        = float(line[38:46]),
                z        = float(line[46:54]),
                occupancy= float(line[54:60]) if line[54:60].strip() else 1.0,
                bfactor  = float(line[60:66]) if line[60:66].strip() else 0.0,
                element  = line[76:78].strip() if len(line) > 77 else "",
                charge   = line[78:80] if len(line) > 79 else "  ",
            )
        except (ValueError, IndexError):
            continue
        # Infer element if missing
        if not a.element:
            a.element = _infer_element(a.name)
        atoms.append(a)
    return atoms


def write_pdb(atoms: List[Atom], path: str) -> None:
    with open(path, "w") as fh:
        for i, a in enumerate(atoms, 1):
            a.serial = i
            fh.write(a.to_pdb_line())
        fh.write("END\n")


def _infer_element(name: str) -> str:
    n = name.strip()
    if not n:
        return ""
    # Strip leading digits and H-disambiguators common in PDB
    for ch in n:
        if ch.isalpha():
            return ch.upper()
    return n[0].upper()

# ---------------------------------------------------------------------------
# Geometry utilities
# ---------------------------------------------------------------------------

Vec3 = Tuple[float, float, float]


def _v(a: Vec3, b: Vec3) -> Vec3:
    """b - a"""
    return (b[0]-a[0], b[1]-a[1], b[2]-a[2])


def _add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])


def _scale(v: Vec3, s: float) -> Vec3:
    return (v[0]*s, v[1]*s, v[2]*s)


def _norm(v: Vec3) -> float:
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def _unit(v: Vec3) -> Vec3:
    n = _norm(v)
    if n < 1e-9:
        return (0.0, 0.0, 0.0)
    return (v[0]/n, v[1]/n, v[2]/n)


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    )


def _dist(a: Vec3, b: Vec3) -> float:
    return _norm(_v(a, b))


def _angle_deg(a: Vec3, b: Vec3, c: Vec3) -> float:
    """Angle at b (degrees)."""
    u = _unit(_v(b, a))
    w = _unit(_v(b, c))
    return math.degrees(math.acos(max(-1.0, min(1.0, _dot(u, w)))))


def rotate_about_axis(p: Vec3, origin: Vec3, axis: Vec3, angle_deg: float) -> Vec3:
    """Rotate point p about axis through origin by angle (degrees)."""
    theta = math.radians(angle_deg)
    u = _unit(axis)
    ux, uy, uz = u
    ct, st = math.cos(theta), math.sin(theta)
    # Translate to origin
    x, y, z = p[0]-origin[0], p[1]-origin[1], p[2]-origin[2]
    # Rodrigues rotation
    dot = ux*x + uy*y + uz*z
    rx = ct*x + st*(uy*z - uz*y) + (1-ct)*dot*ux
    ry = ct*y + st*(uz*x - ux*z) + (1-ct)*dot*uy
    rz = ct*z + st*(ux*y - uy*x) + (1-ct)*dot*uz
    return (rx+origin[0], ry+origin[1], rz+origin[2])


def place_h_tetrahedral_1(center: Vec3, b1: Vec3, b2: Vec3, b3: Vec3,
                           bond_len: float) -> Vec3:
    """
    Place single H in tetrahedral geometry when 3 other substituents are known.
    The H is opposite to the average of the three bond vectors.
    """
    v1 = _unit(_v(center, b1))
    v2 = _unit(_v(center, b2))
    v3 = _unit(_v(center, b3))
    avg = _unit((v1[0]+v2[0]+v3[0], v1[1]+v2[1]+v3[1], v1[2]+v2[2]+v3[2]))
    # Place H opposite to the resultant
    direction = _unit((-avg[0], -avg[1], -avg[2]))
    return _add(center, _scale(direction, bond_len))


def place_h_tetrahedral_2(center: Vec3, b1: Vec3, b2: Vec3, bond_len: float,
                           angle_deg: float = TETRA) -> List[Vec3]:
    """
    Place two H atoms in tetrahedral geometry (methylene).
    b1, b2 are the two known substituents.
    Returns two positions.
    """
    v1 = _unit(_v(center, b1))
    v2 = _unit(_v(center, b2))
    bisect = _unit((v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]))
    # H atoms lie in a plane perpendicular to bisect, at ±half_angle from -bisect
    perp = _unit(_cross(v1, v2))
    if _norm(perp) < 1e-9:
        # v1 and v2 are collinear; pick arbitrary perp
        if abs(v1[0]) < 0.9:
            tmp = (1.0, 0.0, 0.0)
        else:
            tmp = (0.0, 1.0, 0.0)
        perp = _unit(_cross(v1, tmp))

    # Each H is at angle_deg from the bisect-inverted direction
    # Using half the dihedral angle between the two planes
    half_ang = math.radians(angle_deg / 2)
    neg_bisect = (-bisect[0], -bisect[1], -bisect[2])

    pos1 = _add(center, _scale(
        _unit(_add(_scale(neg_bisect, math.cos(half_ang)),
                   _scale(perp, math.sin(half_ang)))), bond_len))
    pos2 = _add(center, _scale(
        _unit(_add(_scale(neg_bisect, math.cos(half_ang)),
                   _scale(perp, -math.sin(half_ang)))), bond_len))
    return [pos1, pos2]


def place_h_tetrahedral_3(center: Vec3, b1: Vec3, bond_len: float,
                           angle_deg: float = TETRA,
                           dihedral_offsets: Tuple[float,...] = (60.0, 180.0, -60.0)
                           ) -> List[Vec3]:
    """
    Place three H atoms (methyl group) around center bonded to b1.
    Returns three positions at staggered dihedrals relative to any reference.
    """
    axis = _unit(_v(b1, center))
    # Build a perpendicular vector
    if abs(axis[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)
    perp = _unit(_cross(axis, ref))

    # The H position template: place at bond_angle from axis, then rotate
    half = math.radians(180.0 - angle_deg)  # angle from -axis
    base_dir = _unit(_add(_scale(axis, math.cos(math.radians(angle_deg))),
                           _scale(perp, math.sin(math.radians(180.0-angle_deg)))))
    # Actually: the H–center–b1 angle = angle_deg
    # direction from center toward H:
    c = math.cos(math.radians(180.0 - angle_deg))
    s = math.sin(math.radians(180.0 - angle_deg))
    base_dir = _add(_scale(axis, c), _scale(perp, s))
    base_dir = _unit(base_dir)

    positions = []
    for offset in dihedral_offsets:
        h_pos = rotate_about_axis(
            _add(center, _scale(base_dir, bond_len)),
            center, axis, offset
        )
        positions.append(h_pos)
    return positions


def place_h_planar(center: Vec3, b1: Vec3, b2: Vec3, bond_len: float) -> Vec3:
    """
    Place H in sp2 plane. center is bonded to b1 and b2 (and H).
    H is in the plane of center,b1,b2, opposite to the bisector of b1–center–b2.
    """
    v1 = _unit(_v(center, b1))
    v2 = _unit(_v(center, b2))
    bisect = _unit((v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]))
    direction = _unit((-bisect[0], -bisect[1], -bisect[2]))
    return _add(center, _scale(direction, bond_len))


def place_h_planar_2(center: Vec3, b1: Vec3, b2: Vec3, bond_len: float,
                     angle_deg: float = PLANAR) -> List[Vec3]:
    """
    Place two H atoms in sp2 geometry (=CH2).
    """
    v1 = _unit(_v(center, b1))
    v2 = _unit(_v(center, b2))
    normal = _unit(_cross(v1, v2))
    bisect = _unit((v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]))

    half = math.radians(angle_deg / 2)
    neg_b = (-bisect[0], -bisect[1], -bisect[2])
    perp_in_plane = _unit(_cross(normal, bisect))

    pos1 = _add(center, _scale(_unit(_add(
        _scale(neg_b, math.cos(half)),
        _scale(perp_in_plane, math.sin(half)))), bond_len))
    pos2 = _add(center, _scale(_unit(_add(
        _scale(neg_b, math.cos(half)),
        _scale(perp_in_plane, -math.sin(half)))), bond_len))
    return [pos1, pos2]


def place_h_linear(center: Vec3, b1: Vec3, bond_len: float) -> Vec3:
    """Place H linearly (sp, e.g., acetylene-like)."""
    direction = _unit(_v(b1, center))
    return _add(center, _scale(direction, bond_len))

# ---------------------------------------------------------------------------
# Residue key mapping
# ---------------------------------------------------------------------------

# Standard backbone atoms
_BB = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}

# For each residue, list of hydrogen placement specs:
#   ("H_name", "center", ["bonded_atom_1", ...], geometry, bond_len, extra_angle)
# geometry: "tet1" | "tet2" | "tet3" | "plan1" | "plan2" | "lin"
#
# Backbone H placement specs (shared, applied per residue context)

def _backbone_specs(resname: str, is_nterm: bool) -> List[dict]:
    """Return backbone hydrogen specs."""
    specs = []
    if is_nterm:
        # NH3+: three H atoms
        specs.append({"names": ["H1","H2","H3"], "center": "N",
                       "bonded": ["CA"], "geo": "tet3", "bl": BL_NH})
    else:
        # Peptide NH — one H, planar
        specs.append({"names": ["H"], "center": "N",
                       "bonded": ["CA", "C_prev"], "geo": "plan1", "bl": BL_NH})

    if resname == "GLY":
        specs.append({"names": ["HA2","HA3"], "center": "CA",
                       "bonded": ["N","C"], "geo": "tet2", "bl": BL_CH_ALI})
    elif resname == "PRO":
        # PRO: no backbone NH; HA only
        specs = [s for s in specs if s["center"] != "N"]
        specs.append({"names": ["HA"], "center": "CA",
                       "bonded": ["N","C","CB"], "geo": "tet1", "bl": BL_CH_ALI})
    else:
        specs.append({"names": ["HA"], "center": "CA",
                       "bonded": ["N","C","CB"], "geo": "tet1", "bl": BL_CH_ALI})
    return specs


# Side-chain hydrogen specs per residue.
# Each spec: dict with keys:
#   names   : list of H atom name(s) to place
#   center  : atom the H is bonded to
#   bonded  : list of atoms bonded to center (used for geometry)
#   geo     : geometry type
#   bl      : bond length (Å)
#   mobile  : True if this group can rotate (OH, NH3, SH, CH3 not in rings)

SC_SPECS: Dict[str, List[dict]] = {

    "ALA": [
        {"names": ["HB1","HB2","HB3"], "center": "CB",
         "bonded": ["CA"], "geo": "tet3", "bl": BL_CH_ALI},
    ],

    "ARG": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG2","HG3"], "center": "CG",
         "bonded": ["CB","CD"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD2","HD3"], "center": "CD",
         "bonded": ["CG","NE"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HE"],  "center": "NE",
         "bonded": ["CD","CZ"], "geo": "plan1", "bl": BL_NH},
        {"names": ["HH11","HH12"], "center": "NH1",
         "bonded": ["CZ"], "geo": "plan2", "bl": BL_NH},
        {"names": ["HH21","HH22"], "center": "NH2",
         "bonded": ["CZ"], "geo": "plan2", "bl": BL_NH},
    ],

    "ASN": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD21","HD22"], "center": "ND2",
         "bonded": ["CG"], "geo": "plan2", "bl": BL_NH,
         "flip_pair": ("OD1","ND2")},
    ],

    "ASP": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
    ],

    "CYS": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","SG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG"],  "center": "SG",
         "bonded": ["CB"], "geo": "lin", "bl": BL_SH,
         "mobile": True},
    ],

    "GLN": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG2","HG3"], "center": "CG",
         "bonded": ["CB","CD"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HE21","HE22"], "center": "NE2",
         "bonded": ["CD"], "geo": "plan2", "bl": BL_NH,
         "flip_pair": ("OE1","NE2")},
    ],

    "GLU": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG2","HG3"], "center": "CG",
         "bonded": ["CB","CD"], "geo": "tet2", "bl": BL_CH_ALI},
    ],

    "GLY": [],  # handled entirely in backbone

    "HIS": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD2"],  "center": "CD2",
         "bonded": ["CG","NE2"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HE1"],  "center": "CE1",
         "bonded": ["ND1","NE2"], "geo": "plan1", "bl": BL_CH_ARO},
        # ND1-H (HID) and NE2-H (HIE) — handled by flip logic
    ],

    "ILE": [
        {"names": ["HB"],  "center": "CB",
         "bonded": ["CA","CG1","CG2"], "geo": "tet1", "bl": BL_CH_ALI},
        {"names": ["HG12","HG13"], "center": "CG1",
         "bonded": ["CB","CD1"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG21","HG22","HG23"], "center": "CG2",
         "bonded": ["CB"], "geo": "tet3", "bl": BL_CH_ALI},
        {"names": ["HD11","HD12","HD13"], "center": "CD1",
         "bonded": ["CG1"], "geo": "tet3", "bl": BL_CH_ALI},
    ],

    "LEU": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG"],  "center": "CG",
         "bonded": ["CB","CD1","CD2"], "geo": "tet1", "bl": BL_CH_ALI},
        {"names": ["HD11","HD12","HD13"], "center": "CD1",
         "bonded": ["CG"], "geo": "tet3", "bl": BL_CH_ALI},
        {"names": ["HD21","HD22","HD23"], "center": "CD2",
         "bonded": ["CG"], "geo": "tet3", "bl": BL_CH_ALI},
    ],

    "LYS": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG2","HG3"], "center": "CG",
         "bonded": ["CB","CD"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD2","HD3"], "center": "CD",
         "bonded": ["CG","CE"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HE2","HE3"], "center": "CE",
         "bonded": ["CD","NZ"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HZ1","HZ2","HZ3"], "center": "NZ",
         "bonded": ["CE"], "geo": "tet3", "bl": BL_NH, "mobile": True},
    ],

    "MET": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG2","HG3"], "center": "CG",
         "bonded": ["CB","SD"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HE1","HE2","HE3"], "center": "CE",
         "bonded": ["SD"], "geo": "tet3", "bl": BL_CH_ALI, "mobile": True},
    ],

    "PHE": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD1"],  "center": "CD1",
         "bonded": ["CG","CE1"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HD2"],  "center": "CD2",
         "bonded": ["CG","CE2"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HE1"],  "center": "CE1",
         "bonded": ["CD1","CZ"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HE2"],  "center": "CE2",
         "bonded": ["CD2","CZ"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HZ"],   "center": "CZ",
         "bonded": ["CE1","CE2"], "geo": "plan1", "bl": BL_CH_ARO},
    ],

    "PRO": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG2","HG3"], "center": "CG",
         "bonded": ["CB","CD"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD2","HD3"], "center": "CD",
         "bonded": ["CG","N"], "geo": "tet2", "bl": BL_CH_ALI},
    ],

    "SER": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","OG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HG"],  "center": "OG",
         "bonded": ["CB"], "geo": "lin", "bl": BL_OH, "mobile": True},
    ],

    "THR": [
        {"names": ["HB"],  "center": "CB",
         "bonded": ["CA","OG1","CG2"], "geo": "tet1", "bl": BL_CH_ALI},
        {"names": ["HG1"], "center": "OG1",
         "bonded": ["CB"], "geo": "lin", "bl": BL_OH, "mobile": True},
        {"names": ["HG21","HG22","HG23"], "center": "CG2",
         "bonded": ["CB"], "geo": "tet3", "bl": BL_CH_ALI},
    ],

    "TRP": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD1"],  "center": "CD1",
         "bonded": ["CG","NE1"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HE1"],  "center": "NE1",
         "bonded": ["CD1","CE2"], "geo": "plan1", "bl": BL_NH},
        {"names": ["HE3"],  "center": "CE3",
         "bonded": ["CD2","CZ3"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HZ2"],  "center": "CZ2",
         "bonded": ["CE2","CH2"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HZ3"],  "center": "CZ3",
         "bonded": ["CE3","CH2"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HH2"],  "center": "CH2",
         "bonded": ["CZ2","CZ3"], "geo": "plan1", "bl": BL_CH_ARO},
    ],

    "TYR": [
        {"names": ["HB2","HB3"], "center": "CB",
         "bonded": ["CA","CG"], "geo": "tet2", "bl": BL_CH_ALI},
        {"names": ["HD1"],  "center": "CD1",
         "bonded": ["CG","CE1"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HD2"],  "center": "CD2",
         "bonded": ["CG","CE2"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HE1"],  "center": "CE1",
         "bonded": ["CD1","CZ"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HE2"],  "center": "CE2",
         "bonded": ["CD2","CZ"], "geo": "plan1", "bl": BL_CH_ARO},
        {"names": ["HH"],   "center": "OH",
         "bonded": ["CZ"], "geo": "lin", "bl": BL_OH, "mobile": True},
    ],

    "VAL": [
        {"names": ["HB"],  "center": "CB",
         "bonded": ["CA","CG1","CG2"], "geo": "tet1", "bl": BL_CH_ALI},
        {"names": ["HG11","HG12","HG13"], "center": "CG1",
         "bonded": ["CB"], "geo": "tet3", "bl": BL_CH_ALI},
        {"names": ["HG21","HG22","HG23"], "center": "CG2",
         "bonded": ["CB"], "geo": "tet3", "bl": BL_CH_ALI},
    ],

    # Nucleic acid residues — simplified backbone ribose/deoxyribose H placement
    "DA": [], "DT": [], "DG": [], "DC": [],
    "A": [],  "U": [],  "G": [],  "C": [],
}

# Histidine protonation states: HIE (NE2-H), HID (ND1-H), HIP (both)
HIS_SPECS_ND1 = {"names": ["HD1"], "center": "ND1",
                  "bonded": ["CG","CE1"], "geo": "plan1", "bl": BL_NH}
HIS_SPECS_NE2 = {"names": ["HE2"], "center": "NE2",
                  "bonded": ["CD2","CE1"], "geo": "plan1", "bl": BL_NH}

# ---------------------------------------------------------------------------
# Residue data structure
# ---------------------------------------------------------------------------

@dataclass
class Residue:
    resname: str
    chain:   str
    resseq:  int
    icode:   str
    atoms:   Dict[str, Atom] = field(default_factory=dict)
    is_nterm: bool = False
    is_cterm: bool = False
    prev_c:  Optional[Vec3] = None   # C of previous residue (for peptide H)


def _build_residues(atoms: List[Atom]) -> List[Residue]:
    """Group atoms into residues in chain order."""
    residues: List[Residue] = []
    current: Optional[Residue] = None
    for a in atoms:
        key = (a.chain, a.resseq, a.icode)
        if current is None or (current.chain, current.resseq, current.icode) != key:
            current = Residue(
                resname=a.resname, chain=a.chain,
                resseq=a.resseq, icode=a.icode
            )
            residues.append(current)
        current.atoms[a.name_stripped] = a
    return residues


def _mark_termini(residues: List[Residue]) -> None:
    """Mark N/C termini per chain."""
    chains: Dict[str, List[int]] = defaultdict(list)
    for i, r in enumerate(residues):
        chains[r.chain].append(i)
    for idxs in chains.values():
        if idxs:
            residues[idxs[0]].is_nterm = True
            residues[idxs[-1]].is_cterm = True
    # Thread prev_c for peptide H placement
    for i in range(1, len(residues)):
        prev = residues[i-1]
        curr = residues[i]
        if prev.chain == curr.chain and "C" in prev.atoms:
            curr.prev_c = prev.atoms["C"].pos

# ---------------------------------------------------------------------------
# Hydrogen placement engine
# ---------------------------------------------------------------------------

def _get_pos(residue: Residue, name: str) -> Optional[Vec3]:
    a = residue.atoms.get(name)
    return a.pos if a else None


def _place_spec(residue: Residue, spec: dict) -> List[Tuple[str, Vec3]]:
    """
    Place one or more H atoms according to a spec dict.
    Returns list of (name, position) tuples.
    """
    center_name = spec["center"]
    bonded_names = spec["bonded"]
    geo = spec["geo"]
    bl = spec["bl"]
    names = spec["names"]

    center = _get_pos(residue, center_name)
    if center is None:
        return []

    # Resolve bonded atom positions; handle "C_prev" for peptide H
    bonded = []
    for bn in bonded_names:
        if bn == "C_prev":
            p = residue.prev_c
        else:
            p = _get_pos(residue, bn)
        if p is not None:
            bonded.append(p)

    results: List[Tuple[str, Vec3]] = []

    if geo == "tet1":
        if len(bonded) < 3:
            return []
        pos = place_h_tetrahedral_1(center, bonded[0], bonded[1], bonded[2], bl)
        results = [(names[0], pos)]

    elif geo == "tet2":
        if len(bonded) < 2:
            return []
        positions = place_h_tetrahedral_2(center, bonded[0], bonded[1], bl, TETRA)
        for i, nm in enumerate(names[:2]):
            if i < len(positions):
                results.append((nm, positions[i]))

    elif geo == "tet3":
        if len(bonded) < 1:
            return []
        positions = place_h_tetrahedral_3(center, bonded[0], bl)
        for i, nm in enumerate(names[:3]):
            if i < len(positions):
                results.append((nm, positions[i]))

    elif geo == "plan1":
        if len(bonded) < 2:
            # Only one bonded atom known — place linearly
            if len(bonded) == 1:
                pos = place_h_linear(center, bonded[0], bl)
                results = [(names[0], pos)]
            return results
        pos = place_h_planar(center, bonded[0], bonded[1], bl)
        results = [(names[0], pos)]

    elif geo == "plan2":
        # Two H on an sp2 N/O with one known bond (NH2 of Asn/Gln/Arg)
        if len(bonded) < 1:
            return []
        # Need a second reference atom to define the plane.
        # For Arg NH1/NH2: bonded[0] = CZ; find NE or another neighbor
        # For Asn ND2: bonded[0] = CG; OD1 is in plane
        # We'll use CG/CZ, look for the carbonyl/guanidinium heavy atom
        ref2 = _find_plane_ref(residue, center_name, bonded[0])
        if ref2 is None:
            # Fall back to placing two H linearly offset
            pos = place_h_linear(center, bonded[0], bl)
            results = [(names[0], pos)]
            return results
        positions = place_h_planar_2(center, bonded[0], ref2, bl, PLANAR)
        for i, nm in enumerate(names[:2]):
            if i < len(positions):
                results.append((nm, positions[i]))

    elif geo == "lin":
        if len(bonded) < 1:
            return []
        pos = place_h_linear(center, bonded[0], bl)
        results = [(names[0], pos)]

    return results


def _find_plane_ref(residue: Residue, center_name: str,
                    bonded1: Vec3) -> Optional[Vec3]:
    """
    For sp2 center, find a second reference atom to define the plane.
    Returns position of second in-plane atom connected through the residue graph.
    """
    center = _get_pos(residue, center_name)
    # Connectivity hints per (center_name) pattern
    hints = {
        "ND2": "OD1",  "NE2": "OE1",
        "NH1": "NE",   "NH2": "NE",
    }
    ref_name = hints.get(center_name)
    if ref_name:
        p = _get_pos(residue, ref_name)
        if p:
            return p
    # Generic: pick any heavy atom in the residue that is close to bonded1
    # (i.e., bonded to the same atom as center) and not the center itself
    best = None
    best_dist = 9999.0
    for a in residue.atoms.values():
        if a.element == "H":
            continue
        if a.name_stripped == center_name:
            continue
        if a.pos == bonded1:
            continue
        d = _dist(a.pos, bonded1)
        if 1.1 < d < 1.7 and d < best_dist:
            best_dist = d
            best = a.pos
    return best


def _is_hydrogen(element: str, name: str) -> bool:
    if element.upper() in ("H", "D"):
        return True
    n = name.strip()
    return n.startswith("H") or n.startswith("D")


def _make_h_atom(name: str, pos: Vec3, template: Atom) -> Atom:
    """Create a new hydrogen Atom from a position and a template heavy atom."""
    pdb_name = _format_h_name(name)
    return Atom(
        record    = template.record,
        serial    = 0,
        name      = pdb_name,
        altloc    = template.altloc,
        resname   = template.resname,
        chain     = template.chain,
        resseq    = template.resseq,
        icode     = template.icode,
        x=pos[0], y=pos[1], z=pos[2],
        occupancy = template.occupancy,
        bfactor   = template.bfactor,
        element   = "H",
        charge    = "  ",
    )


def _format_h_name(name: str) -> str:
    """Format H name into 4-char PDB convention."""
    name = name.strip()
    if len(name) <= 3:
        return f" {name:<3}"
    return f"{name:<4}"

# ---------------------------------------------------------------------------
# Clash / VDW scoring
# ---------------------------------------------------------------------------

def _vdw_radius(element: str) -> float:
    return VDW.get(element.upper(), 1.70)


def _clash_score(h_pos: Vec3, all_atoms: List[Atom],
                 exclude_names: set) -> float:
    """
    Return total clash score for a hydrogen at h_pos.
    Score = sum of (overlap - cutoff) for clashing pairs.
    Negative = fewer clashes.
    """
    r_h = _vdw_radius("H")
    score = 0.0
    for a in all_atoms:
        if a.name_stripped in exclude_names:
            continue
        if a.element.upper() in ("H", "D"):
            continue
        d = _dist(h_pos, a.pos)
        r_sum = r_h + _vdw_radius(a.element)
        overlap = r_sum - d
        if overlap > CLASH_OVERLAP_CUTOFF:
            score += (overlap - CLASH_OVERLAP_CUTOFF) ** 2
    return score


def _hbond_score(h_pos: Vec3, donor_heavy: Vec3, acceptors: List[Atom]) -> float:
    """
    Return negative score for favorable H-bonds (lower = better).
    """
    score = 0.0
    for acc in acceptors:
        da_dist = _dist(donor_heavy, acc.pos)
        if da_dist > HBOND_DIST_CUTOFF:
            continue
        # Angle: donor_heavy — H — acceptor
        v_dh = _unit(_v(donor_heavy, h_pos))
        v_ha = _unit(_v(h_pos, acc.pos))
        cos_a = _dot(v_dh, v_ha)
        angle = math.degrees(math.acos(max(-1.0, min(1.0, cos_a))))
        if angle < HBOND_ANGLE_CUTOFF:
            continue
        # Score by distance (simple Lennard-Jones-like)
        r = _dist(h_pos, acc.pos)
        if r < 1.2:
            continue
        score -= 1.0 / max(r, 0.5)
    return score

# ---------------------------------------------------------------------------
# Mobile hydrogen optimization
# ---------------------------------------------------------------------------

_ACCEPTOR_ELEMENTS = {"O", "N", "S", "F"}


def _optimize_mobile_h(h_name: str, center: Vec3, axis_atom: Vec3,
                        axis: Vec3, init_pos: Vec3,
                        context_atoms: List[Atom],
                        bond_len: float) -> Vec3:
    """
    Rotate a mobile H around the center–axis_atom bond to minimize clashes
    and maximize H-bonds.
    Uses coarse 10° search then fine 1° refinement.
    """
    exclude = set()  # no self-exclusion needed; H not in context yet

    acceptors = [a for a in context_atoms
                 if a.element.upper() in _ACCEPTOR_ELEMENTS
                 and _dist(a.pos, center) < HBOND_DIST_CUTOFF + 1.0]
    heavy_context = [a for a in context_atoms
                     if a.element.upper() not in ("H", "D")]

    best_score = 1e18
    best_angle = 0.0

    # Coarse scan
    for deg in range(0, 360, 10):
        pos = rotate_about_axis(init_pos, center, axis, float(deg))
        s = (_clash_score(pos, heavy_context, exclude)
             + _hbond_score(pos, center, acceptors))
        if s < best_score:
            best_score = s
            best_angle = float(deg)

    # Fine scan around best
    for deg in range(-12, 13):
        angle = best_angle + deg
        pos = rotate_about_axis(init_pos, center, axis, angle)
        s = (_clash_score(pos, heavy_context, exclude)
             + _hbond_score(pos, center, acceptors))
        if s < best_score:
            best_score = s
            best_angle = angle

    return rotate_about_axis(init_pos, center, axis, best_angle)


def _optimize_methyl(positions: List[Vec3], center: Vec3, axis_atom: Vec3,
                     context_atoms: List[Atom]) -> List[Vec3]:
    """
    Optimize methyl group (3 H) by rotating as a rigid body.
    """
    axis = _unit(_v(axis_atom, center))
    heavy_context = [a for a in context_atoms
                     if a.element.upper() not in ("H", "D")]

    best_score = 1e18
    best_offset = 0.0

    for deg in range(0, 120, 6):  # 120° period
        rotated = [rotate_about_axis(p, center, axis, float(deg)) for p in positions]
        score = sum(_clash_score(p, heavy_context, set()) for p in rotated)
        if score < best_score:
            best_score = score
            best_offset = float(deg)

    return [rotate_about_axis(p, center, axis, best_offset) for p in positions]

# ---------------------------------------------------------------------------
# Histidine protonation / flip
# ---------------------------------------------------------------------------

_HIS_RING_ATOMS = {"CG","ND1","CD2","CE1","NE2"}

def _score_his_state(residue: Residue, nd1_h: bool, ne2_h: bool,
                     flipped: bool, context_atoms: List[Atom]) -> float:
    """Score a His protonation + flip state (lower = better)."""
    heavy = [a for a in context_atoms if a.element.upper() not in ("H","D")]
    score = 0.0

    # Penalty for each protonation site that exists
    if nd1_h:
        nd1 = _get_pos(residue, "ND1" if not flipped else "NE2")
        ce1 = _get_pos(residue, "CE1")
        cg  = _get_pos(residue, "CG"  if not flipped else "CD2")
        if nd1 and ce1 and cg:
            h_pos = place_h_planar(nd1, ce1, cg, BL_NH)
            score += _clash_score(h_pos, heavy, set())

    if ne2_h:
        ne2 = _get_pos(residue, "NE2" if not flipped else "ND1")
        cd2 = _get_pos(residue, "CD2" if not flipped else "CE1")
        ce1 = _get_pos(residue, "CE1" if not flipped else "CD2")
        if ne2 and cd2 and ce1:
            h_pos = place_h_planar(ne2, cd2, ce1, BL_NH)
            score += _clash_score(h_pos, heavy, set())

    # Small bonus for HIE (NE2-H) over HID (ND1-H) as default
    if nd1_h and not ne2_h:
        score += 0.05  # slight preference for HIE
    if flipped:
        score += 0.02  # slight preference against flip

    return score


def _determine_his_state(residue: Residue,
                         context_atoms: List[Atom]) -> Tuple[bool, bool, bool]:
    """
    Returns (nd1_protonated, ne2_protonated, flipped).
    Evaluates 6 states: (HID,HIE,HIP) × (unflipped,flipped).
    """
    best = (False, True, False)
    best_score = 1e18

    for nd1_h in (True, False):
        for ne2_h in (True, False):
            if not nd1_h and not ne2_h:
                continue  # must have at least one H
            for flipped in (True, False):
                s = _score_his_state(residue, nd1_h, ne2_h, flipped, context_atoms)
                if s < best_score:
                    best_score = s
                    best = (nd1_h, ne2_h, flipped)

    return best


def _add_his_hydrogens(residue: Residue, context_atoms: List[Atom],
                        optimize: bool) -> List[Atom]:
    """Add protons to histidine ring and optionally optimize state."""
    results = []
    if not optimize:
        # Default: HIE (NE2-H only, unflipped)
        nd1_h, ne2_h, flipped = False, True, False
    else:
        nd1_h, ne2_h, flipped = _determine_his_state(residue, context_atoms)

    # His flipping swaps ND1↔NE2 and CD2↔CE1 atom positions in the ring.
    # For coordinate placement purposes we simply swap which nitrogen gets H.
    # (A true flip would require swapping atom coordinates, which we skip for
    #  simplicity; we only swap H assignment.)

    template_atom = next(iter(residue.atoms.values()))

    if nd1_h:
        nd1 = _get_pos(residue, "ND1")
        ce1 = _get_pos(residue, "CE1")
        cg  = _get_pos(residue, "CG")
        if nd1 and ce1 and cg:
            pos = place_h_planar(nd1, ce1, cg, BL_NH)
            results.append(_make_h_atom("HD1", pos, template_atom))

    if ne2_h:
        ne2 = _get_pos(residue, "NE2")
        cd2 = _get_pos(residue, "CD2")
        ce1 = _get_pos(residue, "CE1")
        if ne2 and cd2 and ce1:
            pos = place_h_planar(ne2, cd2, ce1, BL_NH)
            results.append(_make_h_atom("HE2", pos, template_atom))

    return results

# ---------------------------------------------------------------------------
# Asn / Gln flip optimization
# ---------------------------------------------------------------------------

def _score_flip(residue: Residue, spec: dict, flipped: bool,
                context_atoms: List[Atom]) -> float:
    """Score an Asn/Gln amide orientation (lower = better)."""
    heavy = [a for a in context_atoms if a.element.upper() not in ("H","D")]
    o_name, n_name = spec["flip_pair"]  # e.g. ("OD1","ND2")
    positions = _place_spec(residue, spec)
    if not positions:
        return 0.0
    score = 0.0
    for _, pos in positions:
        score += _clash_score(pos, heavy, set())

    if flipped:
        # Also check if O would clash if roles reversed
        o_pos = _get_pos(residue, o_name)
        n_pos = _get_pos(residue, n_name)
        if o_pos and n_pos:
            # Assess "fake" O position at where N currently is
            score += _clash_score(n_pos, heavy, {n_name}) * 0.5
    return score


def _determine_flip(residue: Residue, spec: dict,
                    context_atoms: List[Atom]) -> bool:
    """Return True if the amide should be flipped for better contacts."""
    s_normal  = _score_flip(residue, spec, False, context_atoms)
    s_flipped = _score_flip(residue, spec, True,  context_atoms)
    # Only flip if meaningfully better (>0.1 score improvement)
    return (s_flipped < s_normal - 0.1)

# ---------------------------------------------------------------------------
# Core hydrogen addition per residue
# ---------------------------------------------------------------------------

def _add_residue_hydrogens(residue: Residue,
                            context_atoms: List[Atom],
                            optimize_mobile: bool,
                            optimize_flips: bool) -> List[Atom]:
    """
    Compute and return hydrogen Atom objects for a single residue.
    context_atoms: nearby atoms from other residues for clash/H-bond scoring.
    """
    resname = residue.resname
    new_h: List[Atom] = []

    # Template atom for new H records
    template = next(iter(residue.atoms.values()))

    def add_h(name: str, pos: Vec3) -> None:
        new_h.append(_make_h_atom(name, pos, template))

    # --- Backbone ---
    is_aa = resname in SC_SPECS or resname in ("GLY",)
    if is_aa and resname != "HOH":
        for spec in _backbone_specs(resname, residue.is_nterm):
            if spec["geo"] == "tet3" and "N" in residue.atoms:
                # N-terminal NH3: treat as mobile (optimize rotation)
                center = _get_pos(residue, "N")
                ca = _get_pos(residue, "CA")
                if center and ca:
                    axis = _unit(_v(ca, center))
                    positions = place_h_tetrahedral_3(center, ca, BL_NH)
                    if optimize_mobile:
                        positions = _optimize_methyl(positions, center, ca, context_atoms)
                    for nm, pos in zip(spec["names"], positions):
                        add_h(nm, pos)
            elif spec["geo"] == "plan1":
                placements = _place_spec(residue, spec)
                for nm, pos in placements:
                    add_h(nm, pos)
            else:
                placements = _place_spec(residue, spec)
                for nm, pos in placements:
                    add_h(nm, pos)

    # --- Side chain ---
    sc_specs = SC_SPECS.get(resname, [])

    for spec in sc_specs:
        is_mobile = spec.get("mobile", False)
        is_flip   = "flip_pair" in spec
        is_methyl = (spec["geo"] == "tet3")

        # Flip optimization for Asn/Gln
        if is_flip and optimize_flips:
            flipped = _determine_flip(residue, spec, context_atoms)
            if flipped:
                # Swap the amide O and N pair — swap names in spec
                o_name, n_name = spec["flip_pair"]
                o_pos = _get_pos(residue, o_name)
                n_pos = _get_pos(residue, n_name)
                if o_pos and n_pos:
                    # Temporarily swap their positions in residue for geometry
                    orig_o = residue.atoms.get(o_name)
                    orig_n = residue.atoms.get(n_name)
                    if orig_o and orig_n:
                        orig_o.x, orig_n.x = orig_n.x, orig_o.x
                        orig_o.y, orig_n.y = orig_n.y, orig_o.y
                        orig_o.z, orig_n.z = orig_n.z, orig_o.z
                        placements = _place_spec(residue, spec)
                        # Restore
                        orig_o.x, orig_n.x = orig_n.x, orig_o.x
                        orig_o.y, orig_n.y = orig_n.y, orig_o.y
                        orig_o.z, orig_n.z = orig_n.z, orig_o.z
                    else:
                        placements = _place_spec(residue, spec)
                else:
                    placements = _place_spec(residue, spec)
            else:
                placements = _place_spec(residue, spec)

            for nm, pos in placements:
                add_h(nm, pos)
            continue

        placements = _place_spec(residue, spec)
        if not placements:
            continue

        if is_methyl and optimize_mobile and len(placements) == 3:
            center = _get_pos(residue, spec["center"])
            bonded_name = spec["bonded"][0]
            axis_atom = _get_pos(residue, bonded_name)
            if center and axis_atom:
                positions = [pos for _, pos in placements]
                positions = _optimize_methyl(positions, center, axis_atom,
                                             context_atoms)
                for nm, pos in zip([n for n, _ in placements], positions):
                    add_h(nm, pos)
                continue

        if is_mobile and len(placements) == 1:
            center = _get_pos(residue, spec["center"])
            bonded_name = spec["bonded"][0]
            axis_atom = _get_pos(residue, bonded_name)
            if center and axis_atom and optimize_mobile:
                h_name, init_pos = placements[0]
                axis = _unit(_v(axis_atom, center))
                opt_pos = _optimize_mobile_h(
                    h_name, center, axis_atom, axis, init_pos,
                    context_atoms, spec["bl"])
                add_h(h_name, opt_pos)
                continue

        for nm, pos in placements:
            add_h(nm, pos)

    # --- Histidine ---
    if resname in ("HIS", "HIE", "HID", "HIP"):
        new_h.extend(_add_his_hydrogens(residue, context_atoms, optimize_flips))

    return new_h

# ---------------------------------------------------------------------------
# Context atom lookup (spatial neighbor search)
# ---------------------------------------------------------------------------

_CONTEXT_RADIUS = 5.0  # Å — atoms within this distance matter for optimization


def _build_context(residue: Residue, all_atoms: List[Atom]) -> List[Atom]:
    """Return atoms within _CONTEXT_RADIUS of any atom in residue, excluding self."""
    own_serials = {a.serial for a in residue.atoms.values()}
    cx, cy, cz = 0.0, 0.0, 0.0
    n = len(residue.atoms)
    if n == 0:
        return []
    for a in residue.atoms.values():
        cx += a.x; cy += a.y; cz += a.z
    cx /= n; cy /= n; cz /= n

    result = []
    for a in all_atoms:
        if a.serial in own_serials:
            continue
        if _is_hydrogen(a.element, a.name):
            continue
        d = math.sqrt((a.x-cx)**2 + (a.y-cy)**2 + (a.z-cz)**2)
        if d < _CONTEXT_RADIUS + 3.0:
            result.append(a)
    return result

# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def add_hydrogens(input_path: str,
                  output_path: str,
                  optimize_mobile: bool = True,
                  optimize_flips:  bool = True,
                  remove_existing: bool = True) -> int:
    """
    Add hydrogens to a PDB structure file.

    Parameters
    ----------
    input_path      : Path to input PDB file.
    output_path     : Path to write output PDB file.
    optimize_mobile : Rotate OH/SH/NH3/methyl to minimize clashes and
                      maximize H-bonds.
    optimize_flips  : Optimize Asn/Gln amide flips and His protonation state.
    remove_existing : Remove any existing hydrogen atoms before adding new ones.

    Returns
    -------
    Number of hydrogen atoms added.
    """
    with open(input_path) as fh:
        lines = fh.readlines()

    # Preserve non-coordinate records
    header_lines = [l for l in lines
                    if not l.startswith(("ATOM  ", "HETATM", "END"))]

    atoms = parse_pdb(lines)

    # Strip existing H atoms
    if remove_existing:
        atoms = [a for a in atoms
                 if not _is_hydrogen(a.element, a.name)]

    residues = _build_residues(atoms)
    _mark_termini(residues)

    # Serially number atoms before context lookup
    for i, a in enumerate(atoms, 1):
        a.serial = i

    added_count = 0
    output_atoms: List[Atom] = list(atoms)

    for residue in residues:
        resname = residue.resname
        # Skip water, ions, and unknown residues without specs
        if resname in ("HOH", "WAT", "H2O"):
            continue
        if resname not in SC_SPECS and resname not in (
                "HIS","HIE","HID","HIP","GLY"):
            continue

        context = _build_context(residue, atoms)

        new_h = _add_residue_hydrogens(
            residue, context, optimize_mobile, optimize_flips
        )
        output_atoms.extend(new_h)
        added_count += len(new_h)

    # Sort: keep original heavy-atom order, insert H after their parent heavy atom
    output_atoms = _sort_atoms(output_atoms, residues)

    # Re-number
    for i, a in enumerate(output_atoms, 1):
        a.serial = i

    with open(output_path, "w") as fh:
        for l in header_lines:
            if l.startswith(("REMARK","TITLE","CRYST","SCALE","ORIGX","SEQRES")):
                fh.write(l)
        for a in output_atoms:
            fh.write(a.to_pdb_line())
        fh.write("END\n")

    return added_count


def _sort_atoms(atoms: List[Atom], residues: List[Residue]) -> List[Atom]:
    """
    Sort atoms so that within each residue, H atoms appear after
    the heavy atom they are bonded to (heuristic: by name prefix).
    """
    # Build residue ordering key: (chain, resseq, icode)
    res_order: Dict[Tuple, int] = {}
    for i, r in enumerate(residues):
        key = (r.chain, r.resseq, r.icode)
        res_order[key] = i

    def atom_key(a: Atom):
        rkey = (a.chain, a.resseq, a.icode)
        ri = res_order.get(rkey, 9999)
        is_h = 1 if _is_hydrogen(a.element, a.name) else 0
        return (ri, is_h, a.name_stripped)

    return sorted(atoms, key=atom_key)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main() -> None:
    parser = argparse.ArgumentParser(
        description="Add hydrogens to a protein PDB structure file. "
                    "Reimplementation of the Richardson lab 'reduce' utility."
    )
    parser.add_argument("input",  help="Input PDB file")
    parser.add_argument("output", help="Output PDB file")
    parser.add_argument("--no-optimize-mobile", dest="optimize_mobile",
                        action="store_false", default=True,
                        help="Skip rotation optimization of OH/SH/NH3/methyl groups")
    parser.add_argument("--no-flip", dest="optimize_flips",
                        action="store_false", default=True,
                        help="Skip Asn/Gln flip and His protonation optimization")
    parser.add_argument("--keep-hydrogens", dest="remove_existing",
                        action="store_false", default=True,
                        help="Keep existing hydrogens (do not strip before adding)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading {args.input} ...")
    n = add_hydrogens(
        args.input, args.output,
        optimize_mobile  = args.optimize_mobile,
        optimize_flips   = args.optimize_flips,
        remove_existing  = args.remove_existing,
    )
    print(f"Added {n} hydrogen atoms -> {args.output}")


if __name__ == "__main__":
    _main()

#!/usr/bin/env python
# Parses output of FM and regurgitates in the form of a gromacs itp file.

import sys, os
dn = os.path.dirname
sys.path.append(dn(dn(os.path.abspath(__file__))))

from mol import read_mol
from cg_topol import *
from numpy import array, sqrt, sum, argmax, newaxis, abs, dot, pi, arange
from read_top import TOP

def main(argv):
    assert len(argv) == 4, "Usage: %s <mol.sdf> <param dir> <out.itp>"%argv[0]
    mol = read_mol(argv[1])
    top = top_of_param(argv[2], mol)
    top.write(argv[3])

def get_term(name):
    id, c = read_poly_term(name)
    return tuple(map(tname, id[1].split("-"))), c

# Insert optional name mangling scheme here:
def tname(n):
    return n

def top_of_param(path, mol):
    # fn type, comb. rule, gen-pairs, fudgeLJ, fudgeCoul
    defaults = { 'default':[1, 2, True, 1.0, 0.75] }
    atoms = {}
    for n in sorted(set(mol.t)):
        a = mol.t.index(n)
        atoms[tname(n)] = 1, mol.m[a], mol.q[a], 'A', 0.1, 0.0

    bonds = {}
    constraints = {}
    angles = {}
    pair14 = {}
    dihedrals = {}
    nonbonded = {}
    to_deg = 180./pi
    for fname in os.listdir(path):
	tp = fname.split("_")[0]
	name = os.path.join(path, fname)
	if tp == "pbond":
	    id, c = get_term(name)
	    bonds[id] = (1, -0.5*c[1]/c[2], 2*c[2])
	elif tp == "pangle":
	    id, c = get_term(name)
	    if angles.has_key(id): # handle case where UB is present
		r = angles[id][3:]
	    else:
		r = ()
	    angles[id] = (1+4*(len(r) > 0), -to_deg*0.5*c[1]/c[2], 2*c[2]) + r
	elif tp == "pub":
	    id, c = get_term(name)
	    if angles.has_key(id):
                r = angles[id][1:]
	    else:
		r = 100.0, 0.0
	    angles[id] = (5, r[0], r[1], -0.5*c[1]/c[2], 2*c[2])
	elif tp == "ptor":
	    id, c = get_term(name)
            # connecting to MMFF94,
            # V3 = c[3]/4.0
            # V2 = -0.5*c[2]
            # V1 = c[1] + 3*c[3]/4.0
            # V0 = c[1] - c[2] + c[3]
            const = c[0] + c[1] - c[2] + c[3]
	    dihedrals[id] = (3, const, -c[1], c[2], -c[3], c[4], 0.0)
	elif tp == "pimprop":
	    # cg_topol/pimprop.py:117
	    # just writes the line, "#IMPR <name> <K>"
	    line = open(name).read()[0].split()
	    id = map(tname, line[1].split("_")[1].split("-"))
	    dihedrals[tuple(id)] = (2, 0.0, float(line[2]))
	elif tp == "ljpair":
	    # pair_terms name = "4+%s-%s"
	    id, c = get_term(name)
	    c6, c12 = c[1:]
	    eps, R0 = 0.25*c6*c6/c12, (-2*c12/c6)**(1/6.0)
	    if "4+" not in id[0]:
                raise ValueError, "Invalid LJ format."
            num = id[0].split('+')
            id = (tname(num[1]),) + id[1:]
            nonbonded[id] = (R0, eps)
    return TOP(defaults, atoms, bonds, constraints, angles,
               dihedrals, pair14, nonbonded)

if __name__=="__main__":
    main(sys.argv)


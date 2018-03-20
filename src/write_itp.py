#!/usr/bin/env python
# Parses output of FM and regurgitates in the form of a gromacs itp file.

import sys, os
from mol import read_mol
from psf import read_psf
from forcesolve import *
from numpy import array, sqrt, sum, argmax, abs, dot, pi, arange
from top import TOP

def main(argv):
    assert len(argv) == 3, "Usage: %s <param dir> <out.itp>"%argv[0]
    top = top_of_param(argv[1])
    top.write(argv[2])

# necessary scales to convert to nm, kJ/mol
dist = 0.1
en   = 4.184

def get_term(name):
    id, c = read_poly_term(name)
    return tuple(map(tname, id[1:])), c

# Insert optional name mangling scheme here:
def tname(n):
    return n

def top_of_param(path):
    mol = read_psf(os.path.join(path, "molecule.psf"))
    # fn type, comb. rule, gen-pairs, fudgeLJ, fudgeCoul
    #defaults = { 'default':[1, 2, True, 1.0, 0.75] }
    defaults = {}
    atoms = {}
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
	    bonds[id] = (1, -0.5*c[1]/c[2]*dist, 2*c[2]*en/dist**2)
	elif tp == "pangle":
	    id, c = get_term(name)
	    if angles.has_key(id): # handle case where UB is present
		r = angles[id][3:]
	    else:
		r = ()
	    angles[id] = (1+4*(len(r) > 0), -to_deg*0.5*c[1]/c[2], 2*c[2]*en) \
                       + r
	elif tp == "pub":
	    id, c = get_term(name)
	    if angles.has_key(id):
                r = angles[id][1:]
	    else:
		r = 100.0, 0.0
            r0 = -0.5*c[1]/c[2]*dist
            k  = 2*c[2]*en/dist**2
            if k > 1e-6: # was zeroed during fitting
                angles[id] = (5, r[0], r[1], r0, k)
	elif tp == "ptor":
	    id, c = get_term(name)
            c *= en
            # connecting to MMFF94,
            # V3 = c[3]/4.0
            # V2 = -0.5*c[2]
            # V1 = c[1] + 3*c[3]/4.0
            # V0 = c[1] - c[2] + c[3]
            for i in range(1, len(c), 2): # switch zero to 180 deg.
                c[i] *= -1.0
            c[0] = -sum(c[1:]) # set E(0 deg at x = 1) = 0
            dihedrals[id + (3,)] = tuple(c)
	elif tp == "pimprop":
	    # cg_topol/pimprop.py:117
	    # just writes the line, "#IMPR <name> <K>"
	    line = open(name).readlines()[0].split()
            id = map(tname, line[1].split("_")[1:])
	    dihedrals[tuple(id) + (2,)] = (0.0, float(line[2])*en)
	elif tp == "ljpair":
	    # pair_terms name = "4+%s-%s"
	    id, c = get_term(name)
	    c6, c12 = c[1:]
	    eps, R0 = 0.25*c6*c6/c12, (-2*c12/c6)**(1/6.0)
	    if "4+" not in id[0]:
                raise ValueError, "Invalid LJ format."
            num = id[0].split('+')
            id = (tname(num[1]),) + id[1:]
            if eps < 1e-9:
                eps = 0.0
                R0 = 3.0
            nonbonded[id] = (R0*dist, eps*en)

    for n in sorted(set(mol.t)):
        a = mol.t.index(n)
        try:
            nb = nonbonded[(n,n)]
        except KeyError:
            nb = (0.3, 0.0)
        atoms[tname(n)] = (1, mol.m[a], mol.q[a], 'A') + nb

    return TOP(defaults, atoms, bonds, constraints, angles,
               dihedrals, pair14, nonbonded)

if __name__=="__main__":
    main(sys.argv)


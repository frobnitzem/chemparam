#!/usr/bin/env python2.7
# Driver program for processing Gromacs/CHARMM-specific
# information into PDB and cg_topol data structures
# used by forcesolve.
#
# This program assumes units of:
# en ~ kcal/mol
# dist ~ Angstroms
# charge ~ e_0

import sys, os, argparse
from numpy import load, newaxis, sum, zeros, newaxis
from mol import read_mol, onefour
from top import read_top
from psf import read_psf
from itp import read_itp
from prm import read_prm
from forcesolve import *

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Create Molecule Topology')
    parser.add_argument('mol', metavar='sys.mol', type=str,
			help='System MOL/PSF file.')
    parser.add_argument('out', metavar='out.itp', type=str,
			help='Output itp name.')
    # Toggles.
    parser.add_argument('--noUB', dest='UB', action='store_false',
		        help='Don\'t Use Urey-Bradley Terms.')
    parser.add_argument('--14', dest='pairs', action='store_true',
		        help='Fit special 14 pairs.')
    return parser.parse_args()

def main(args):
    out = args.out
    # Read input data.
    if args.mol[-4:] == '.psf':
        mol = read_psf(args.mol)
    else:
        mol = read_mol(args.mol)
    mol.write_itp(out, os.path.basename(args.mol).rsplit('.',1)[0])

if __name__=="__main__":
    args = parse_args(sys.argv)
    print(args)
    main(args)


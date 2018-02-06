from mol import *
import sys

if __name__=="__main__":
    assert len(sys.argv) == 4, \
           "Usage: %s <in.mol> <out.psf> <out.itp>"%sys.argv[0]
    m = read_mol(sys.argv[1])
    m.write_psf(sys.argv[2])
    m.write_itp(sys.argv[3])


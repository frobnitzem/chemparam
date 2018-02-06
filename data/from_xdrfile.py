#!/usr/bin/env python2.7

import sys
from numpy import save, array
from xdrfile import xdrfile

def main(argv):
    assert len(argv) == 3, "Usage: %s <in.trr> <out.xf.npy>"%argv[0]
    t = xdrfile(argv[1])
    xf = []
    for f in t:
        xf.append([f.x.copy(), f.f.copy()])

    # convert to Ang, kcal/mol/Ang
    xf = array(xf)
    xf[:,0] *= 10.0
    xf[:,1] *= 0.1/4.184
    save(argv[2], xf)
    print("Remember to check for artifacts from periodicity"
          " (not corrected here).")

if __name__ == "__main__":
    main(sys.argv)


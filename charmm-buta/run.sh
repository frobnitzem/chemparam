echo "Simple match of butane parameters in CHARMM-style FF."

rm -fr out
mkdir out
../src/list_terms.py buta.psf buta.itp
../src/match.py --param cgenff.36.prm buta.itp mm-buta.npy out
../src/write_prm.py out out/out.prm

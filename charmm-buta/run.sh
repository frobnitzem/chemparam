echo "Simple match of butane parameters in CHARMM-style FF."

rm -fr out
mkdir out
../src/match.py --top cgenff.36.prm buta.psf mm-buta.npy out
../src/write_prm.py out out/out.prm

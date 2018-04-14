echo "Simple match of butane parameters in CHARMM-style FF."

rm -fr out
mkdir out
../src/match.py ../data/charmm_types.py buta.psf mm-buta.npy out
#../src/match.py --param ../data/cgenff.36.prm --noLJ charmm-noLJ.py buta.psf mm-buta.npy out
../src/write_prm.py out out/out.prm

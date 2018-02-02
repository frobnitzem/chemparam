#!/bin/sh

if [ $# -ne 4 ]; then
    echo "Usage: $0 <in.mol> <nconf> <theory> <out.itp>"
    exit 1
fi

rm -fr out perm scratch
mkdir -p out perm scratch

DIR=$PWD

name=`basename $1 .mol`
# 1. Generate force data for 1000 conformers.
$DIR/gen_frc.py $1 $2 $3 $name.xf.npy
# 2. Perform FM with no charge or LJ fitting.
$DIR/match.py --top $DIR/nonbond.itp --noLJ $1 $name.xf.npy out
# 3. Write the topology in Gromacs format.
$DIR/write_top.py $1 out $4

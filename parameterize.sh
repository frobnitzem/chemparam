#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: $0 <in.mol> <nconf> <theory>"
    exit 1
fi

rm -fr out perm scratch
mkdir -p out perm scratch

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

name=`basename $1 .mol`
# 1. Generate force data.
$DIR/src/gen_frc.py $1 $2 $3 $name.xf.npy
# 2. Perform FM with no charge or LJ fitting.
$DIR/src/match.py --top $DIR/nonbond.itp --noLJ $1 $name.xf.npy out
# 3. Write the topology in Gromacs format.
$DIR/src/write_top.py $1 out $name.itp
# 4. Write the topology in Charmm format.
$DIR/src/write_prm.py out $name.prm

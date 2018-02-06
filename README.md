# ChemParam

chemparam - an automated interface for calling [forcesolve](https://github.com/frobnitzem/forcesolve)
using Gromacs and CHARMM-style interaction energy forms.

Although it's been available since 2008 as the first general-purpose
force matching software ever, few research groups have developed methods
which use it.  Forcesolve is blindingly fast, rock solid, and can
fit arbitrary generalized additive forms.

This code uses OpenBabel to read in molecules and assign
initial chemical atom types, sample molecular conformers,
then (optionally) libnwchem to compute forces.  Next,
it runs forcesolve to match parameters to the best of the
fitting form's ability and outputs them in both Gromacs `itp`
and CHARMM `prm` formats.

## Controllable Options

1. Can fit pairwise LJ terms or assign from an existing parameter file.
2. Can either fit charges or leave them alone.
3. Can include Urey-Bradley terms.
4. Can fit special 1-4 terms or assign the same LJ terms everywhere.
5. Can handle special 1-4 LJ or charge scaling.
   * Caveat: Cannot fit LJ using a special 1-4 scale.

## Dependencies

* [Python 2 or 3](https://www.python.org)
* [OpenBabel 2.3](https://openbabel.org)
** requires pybel option: `apt-get install python-openbabel`
* [ForceSolve](https://github.com/frobnitzem/forcesolve)
** Copy it to any site-packages, there's no C to compile there
* Optional: [LibNWChem](https://github.com/nwchemgit/nwchem/pull/13)
** Use NWChem to compute your forces.
** If you write a wrapper to [Psi4](http://www.psicode.org/) or [PySCF](http://sunqm.github.io/pyscf), let me know!

## Running

The `parameterize.sh` script should be a one-shot wonder.  Try it out on your favorite `mol` file while enjoying as many handfuls of popcorn as possible (but you probably won't get past 2).

Most of the work is done by `match.py`, which has self-explanatory long-style commandline options.  No extra assembly required, no BS.

## Contributing

Pull requests and bug reports are welcome!

## Authors

* **David M. Rogers** [PredictiveStatMech](https://predictivestatmech.org)

## Acknowledgments

* This code was made possible by the persistence of Phillip S. Hudson and H. Lee Woodcock.

# IF YOU USE THIS SOFTWARE

Please cite ChemParam as:
> David M. Rogers, ChemParam (2018) https://github.com/frobnitzem/chemparam.


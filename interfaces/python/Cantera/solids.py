"""This module defines classes and functions used to model gas mixtures."""

from solution import Solution

def Solid(src="",
                kmodel=1, transport=None):
    return Solution(import_file=import_file,
                    thermo_db="",
                    eos=0,
                    id=id,
                    kmodel=kmodel,
                    trmodel=transport,
                    validate=0)

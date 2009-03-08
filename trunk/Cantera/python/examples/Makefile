#!/bin/sh

PY_DEMO_DIRS = equilibrium flames gasdynamics kinetics liquid_vapor \
               misc reactors surface_chemistry transport fuel_cells

run-all:
	@(for d in $(PY_DEMO_DIRS) ; do \
	    echo "entering directory $${d}..."; \
            (cd $${d}; make run) ; \
	done)

cleanup:
	@(for dd in $(PY_DEMO_DIRS) ; do \
	    echo "entering directory $${dd}..."; \
            (cd $${dd}; make -i clean; cd ..) ; \
	done)

# end of file

# 
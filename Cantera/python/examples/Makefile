#!/bin/sh

PY_DEMO_DIRS = equilibrium flames gasdynamics kinetics liquid_vapor \
               misc reactors surface_chemistry transport fuel_cells
MAKE = gmake
all:
	@(for d in $(PY_DEMO_DIRS) ; do \
            echo "entering directory $${d}..."; \
            (cd $${d}; $(MAKE)) ; \
        done)

run:
	@(for d in $(PY_DEMO_DIRS) ; do \
	    echo "entering directory $${d}..."; \
            (cd $${d}; $(MAKE) run) ; \
	done)

test:
	@(for d in $(PY_DEMO_DIRS) ; do \
            echo "entering directory $${d}..."; \
            (cd $${d}; $(MAKE) test) ; \
        done)

install:
	/usr/bin/install -c -d  /home/hkmoffa/arch/linux/cantera-1.7_develop/demos/python
	/usr/bin/install -c -c -m ug+rw,o+r Makefile /home/hkmoffa/arch/linux/cantera-1.7_develop/demos/python 
	/usr/bin/install -c -c -m ug+rw,o+r run_examples.py /home/hkmoffa/arch/linux/cantera-1.7_develop/demos/python 
	@(for d in $(PY_DEMO_DIRS) ; do \
            echo "entering directory $${d}..."; \
            (cd $${d}; $(MAKE) install) ; \
        done)


clean:
	@(for dd in $(PY_DEMO_DIRS) ; do \
	    echo "entering directory $${dd}..."; \
            (cd $${dd}; $(MAKE) -i clean; cd ..) ; \
	done)

# end of file



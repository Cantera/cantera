#!/bin/sh

PY_DEMOS = flame1.py flame2.py stflame1.py npflame1.py free_h2_air.py \
           adiabatic_flame.py fixed_T_flame.py

run:
	@(for py in $(PY_DEMOS) ; do \
	    echo "running $${py}..."; \
          $(PYTHON_CMD) "$${py}"; \
	done)

clean:
	rm -f *.log *.csv *.xml

# end of file


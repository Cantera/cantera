#!/bin/sh

PY_DEMOS = critProperties.py rankine.py
PYTHON_CMD = /usr/bin/python

run:
	@(for py in $(PY_DEMOS) ; do \
	    echo "running $${py}..."; \
          $(PYTHON_CMD) "$${py}"; \
	done)

clean:
	rm -f *.log *.csv *.xml

# end of file


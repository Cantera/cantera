#!/bin/sh

PY_DEMOS = dustygas.py
PYTHON_CMD = /Library/Frameworks/Python.framework/Versions/Current/bin/python

run:
	@(for py in $(PY_DEMOS) ; do \
	    echo "running $${py}..."; \
          $(PYTHON_CMD) "$${py}"; \
	done)

clean:
	rm -f *.log *.csv *.xml

# end of file


#!/bin/sh

PY_DEMOS = ratecoeffs.py

run:
	@(for py in $(PY_DEMOS) ; do \
	    echo "running $${py}..."; \
          $(PYTHON_CMD) "$${py}"; \
	done)

clean:
	rm -f *.log *.csv *.xml

# end of file


#!/bin/sh

PY_DEMOS = combustor.py function1.py mix1.py mix2.py piston.py reactor1.py \
           reactor2.py sensitivity1.py

run:
	@(for py in $(PY_DEMOS) ; do \
	    echo "running $${py}..."; \
          $(PYTHON_CMD) "$${py}"; \
	done)

clean:
	rm -f *.log *.csv *.xml

# end of file


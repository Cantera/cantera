import sys
import string
from distutils.core import setup, Extension

import ctconf

platform = sys.platform

# values:
#  0   do nothing
#  1   install only ctml_writer.py
#  2   install full package

if ctconf.buildPython >= 2:
    setup(name="Cantera",
          version=ctconf.ctversion,
          description="The Cantera Python Interface",
          long_description="""
              """,
          author="Prof. D. G. Goodwin, Caltech",
          author_email="dgoodwin@caltech.edu",
          url="http://code.google.com/p/cantera",
          package_dir = {'MixMaster':'../../apps/MixMaster'},
          packages = ["","Cantera","Cantera.OneD",
                      "MixMaster","MixMaster.Units"],
          package_data = {'Cantera': ['_cantera.so']})
elif ctconf.buildPython == 1:
    setup(name="Cantera CTI File Processor",
          version=ctconf.ctversion,
          description="Converts .cti files to CTML",
          long_description="""
              """,
          author="Prof. D. G. Goodwin, Caltech",
          author_email="dgoodwin@caltech.edu",
          url="http://www.cantera.org",
          py_modules = ["ctml_writer"])

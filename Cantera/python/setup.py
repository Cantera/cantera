#!/usr/bin/python2

from distutils.core import setup, Extension

libs = []
import sys
platform = sys.platform

if platform == "win32":
    libs = ["cantera14"]
else:
    libs = ["ct14", "zeroD","oneD","converters", "transport",
            "cantera","recipes","ctlapack",
            "ctblas", "ctmath", "cvode", "stdc++", "g2c", "m"]


#if sys.argv[1] == 'install':
#    import time
#    f = open('_pydate.py','w')
#    f.write('date = '+`time.time()`)
#    f.close()

setup(name="Cantera",
      version="1.4",
      description="The Cantera Python Interface",
      long_description="""
      """,
      author="Prof. D. G. Goodwin, Caltech",
      author_email="dgoodwin@caltech.edu",
      url="http://www.cantera.org",
      package_dir = {'MixMaster':'../../apps/MixMaster'},
      packages = ["","Cantera","MixMaster","MixMaster.Units"],
      ext_modules=[
          Extension("Cantera._cantera",
                    ["src/pycantera.cpp"],
                    include_dirs=["../../build/include",
                                  "src", "../clib/src"],
                    library_dirs = ["../../build/lib"], libraries = libs)
          ],
     )


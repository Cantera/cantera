from distutils.core import setup, Extension

libdir = ['../../build/lib/i686-pc-win32']
    
setup(name="Cantera",
          version="1.7.0",
          description="The Cantera Python Interface",
          long_description="""
          """,
          author="Prof. D. G. Goodwin, Caltech",
          author_email="dgoodwin@caltech.edu",
          url="http://www.cantera.org",
          package_dir = {'MixMaster':'../../apps/MixMaster'},
          py_modules = ["ctml_writer"],
          packages = ["Cantera","Cantera.OneD",
                      "MixMaster","MixMaster.Units"],
          ext_modules=[ Extension("Cantera._cantera",
                                  ["src/pycantera.cpp"],
                                  include_dirs=["../../build/include",
                                                "src", "../clib/src"],
                                  library_dirs = libdir,
                                  depends = [libdir[0]+'/clib.lib'],
                                  libraries = ["clib"])
                        ],
          )


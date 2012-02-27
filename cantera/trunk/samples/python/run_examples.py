import os
from os.path import walk
pycmd = os.getenv("PYTHON_CMD")
if not pycmd:
    pycmd = "python"

def run_example(dir, file):
    print "******************************************"
    print "  Example "+file+" ("+dir+")"
    print "******************************************"
    r = os.system("cd "+dir+"; "+pycmd+" "+file)
    if r <> 0:
        print "error!"

def run_examples(a, dir, files):
    for f in files:
        base, ext = os.path.splitext(f)
        print base, " <> ",ext
        if dir <> "." and ext == ".py":
            run_example(dir, f)

walk(".",run_examples,None)

## run_example("equilibrium","adiabatic.py")
## run_example("equilibrium","multiphase_plasma.py")
## run_example("equilibrium","plotting.py")
## run_example("equilibrium","simple.py")
## run_example("equilibrium","stoich.py")
## run_example("flames","adiabatic_flame.py")
## run_example("flames","flame1.py")
## run_example("flames","flame2.py")
## run_example("flames","free_h2_air.py")
## run_example("flames","npflame1.py")
## run_example("flames","stflame1.py")
## run_example("gasdynamics","isentropic.py")
## run_example("kinetics","ratecoeffs.py")
## run_example("liquid_vapor","critProperties.py")
## run_example("liquid_vapor","rankine.py")
## run_example("misc","rxnpath1.py")
## run_example("reactors","function1.py")
## run_example("reactors","mix1.py")
## run_example("reactors","mix2.py")
## run_example("reactors","reactor1.py")
## run_example("reactors","reactor2.py")
## run_example("surface_chemistry","catcomb.py")
## run_example("surface_chemistry","diamond.py")
## run_example("transport","dustygas.py")

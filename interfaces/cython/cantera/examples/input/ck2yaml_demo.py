"""
Illustrate conversion of chemkin input files to cantera yaml input files.

The script runs the conversion utility `ck2yaml` as a shell command. The
following GRI-3.0 derived input files are re-created:
 * gri30.yaml (without deprecated phases)
 * h2o2.yaml
 * air.yaml
 * argon.yaml
 * airNASA9.yaml
 * nasa.yaml

Requires: cantera >= 2.5.0
"""

import subprocess

opt = {}

# options for gri30.yaml
opt['gri30'] = ["--input=gri30.inp", "--output=gri30.yaml",
                "--thermo=gri30_thermo.dat", "--transport=gri30_tran.dat",
                "--name=gri30", "--bibtex=gri30.bib"]

# options for h2o2.yaml
opt['h2o2'] = ["--input=h2o2.inp", "--output=h2o2.yaml",
               "--transport=gri30_tran.dat",
               "--name=ohmech"]

# options for air.yaml
opt['air'] = ["--input=air.inp", "--output=air.yaml",
              "--transport=gri30_tran.dat",
              "--name=air"]

# options for argon.yaml (the --quiet option suppresses warnings about
# unexpected/unused species in thermo data)
opt['argon'] = ["--input=argon.inp", "--output=argon.yaml",
                "--thermo=gri30_thermo.dat", "--transport=gri30_tran.dat",
                "--name=argon", "--quiet"]

# options for airNASA9.yaml
opt['airNASA9'] = ["--input=airNASA9.inp", "--output=airNASA9.yaml",
                   "--thermo=airDataNASA9.dat",
                   "--name=airNASA9", "--bibtex=nasa9.bib"]

# options for nasa.yaml
opt['nasa'] = ["--thermo=nasathermo.dat", "--output=nasa.yaml",
               "--bibtex=nasa7.bib", "--permissive"]

# display help text
print('Executing shell command:')
print('ck2yaml --help\n')
subprocess.run(['ck2yaml', '--help'])

for o in opt:

    # print command
    print('\nExecuting shell command to create `{}.yaml`:'.format(o))
    print('ck2yaml {}\n'.format(' '.join(opt[o])))
    subprocess.run(['ck2yaml'] + opt[o])

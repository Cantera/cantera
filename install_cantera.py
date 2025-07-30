# Small python script for Cantera CERFACS installation

import sys
import os
import subprocess


def fill_cantera_conf(file_path, install_dir_path):
    f = open(file_path, 'w')

    text = """\
prefix = '{0}'
"""

    f.write(text.format(install_dir_path))

    if argument == 'local':
        text = """\
boost_inc_dir = '/opt/homebrew/Cellar/boost/{0}/include'
verbose_tests = True
"""
        # use the output of ls /opt/homebrew/Cellar/boost/ to get the version number
        f.write(text.format(subprocess.check_output('ls /opt/homebrew/Cellar/boost/', shell=True).decode('utf-8').split('\n')[0]))
    else :
        if cluster == 'CALYPSO':
            text = """\
boost_inc_dir = '/softs/local/boost/1.86.0/include'
"""
        elif cluster == 'KRAKEN':
            text = """\
boost_inc_dir = '/softs/local/boost/1.78.0_gcc112/include'
"""
        f.write(text)

    f.close()


def execute_with_live_display(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)

    while True:
        # Read a line from the output
        line = process.stdout.readline()

        # Break the loop if there are no more lines
        if not line and process.poll() is not None:
            break

        # Display the line in real-time
        print(line.strip())

    # Wait for the process to complete
    process.wait()

    # Return the exit code of the process
    return process.returncode


# Asking about the environment
if sys.version_info[0] < 3:
    quit('You should compile cantera with python 3 ! Use module load python/3.9.5')
    string_argument = raw_input('Are you installing Cantera on NFS machine ? (yes/no) ')
    python_version = 2
else:
    string_argument = input('Are you installing Cantera on NFS machine ? (yes/no) ')
    python_version = 3

if string_argument in ['yes', 'y']:
    argument = 'NFS'
    cluster = input('On which cluster are you installing Cantera ? (KRAKEN/CALYPSO) ')
elif string_argument in ['no', 'n']:
    argument = 'local'
else:
    quit('Invalid answer ! \n Really ?! How could you fail a yes or no question ?')

dir_path = os.path.dirname(os.path.realpath(__file__))

install_dir_path = os.path.dirname(os.path.realpath(__file__)) + "/INSTALL_DIR"

if not os.path.isdir(install_dir_path):
    os.mkdir(install_dir_path)

#create mech_lib folder in INSTALL_DIR
if not os.path.isdir(install_dir_path + "/mech_lib"):
    os.mkdir(install_dir_path + "/mech_lib")

file_path = dir_path + "/cantera.conf"

if argument == 'local':
    print("""To run cantera, the installation of: \n
    - python, boost and gcc with brew command. \n
    - numpy, cython, scons, wheel and ruamel.yaml with pip3 command. \n
    are required.""")
    update_argument = input('Do you want to install/update those libraries ? (yes/no) ')

    if update_argument in ['yes', 'y']:
        if subprocess.call('which brew', shell=True) == 1:
            subprocess.call('/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"', shell=True)
        subprocess.call('brew install boost', shell=True)
        subprocess.call('brew install gcc', shell=True)
        #subprocess.call('pip3 install numpy --no-use-pep517', shell=True)
        subprocess.call('pip3 install numpy', shell=True)
        subprocess.call('pip3 install cython==3.1.1', shell=True)
        subprocess.call('pip3 install packaging', shell=True)
        # subprocess.call('pip3 install scons==3.1.2', shell=True)
        subprocess.call('pip3 install scons', shell=True)
        subprocess.call('pip3 install wheel', shell=True)
        subprocess.call('pip3 install ruamel.yaml', shell=True)
    elif update_argument in ['no', 'n']:
        pass
    else:
        quit('Invalid answer ! \n Really ?! How could you fail a yes or no question ?')

    fill_cantera_conf(file_path, install_dir_path)
    error = subprocess.call('scons build && scons install', shell=True)
else:
    print("""To run cantera, the installation of: \n
    - numpy, cython, scons, wheel and ruamel.yaml with pip command. \n
    are required.""")
    update_argument = input('Do you want to install/update those libraries ? (yes/no) ')

    if update_argument in ['yes', 'y']:
        subprocess.call('pip3 install numpy', shell=True)
        subprocess.call('pip3 install cython==3.1.1', shell=True)
        subprocess.call('pip3 install packaging', shell=True)
        subprocess.call('pip3 install scons', shell=True)
        subprocess.call('pip3 install wheel', shell=True)
        subprocess.call('pip3 install ruamel.yaml', shell=True)
    elif update_argument in ['no', 'n']:
        pass
    else:
        quit('Invalid answer ! \n Really ?! How could you fail a yes or no question ?')


    fill_cantera_conf(file_path, install_dir_path)
    error = execute_with_live_display("./run_compil_" + cluster)


##only print the following if the installation was successful
#if argument == 'NFS' and not error:
    ##add an alias to load necessary modules
    #print("We suggest to add the alias 'load_cantera_mod' to purge and load necessary modules in ~/.bashrc")
    ##ask user if he wants to add the alias to ~/.bashrc
    #add_alias = input('Do you want to add the alias to ~/.bashrc ? (yes/no) ')
    # if add_alias in ['yes', 'y']:
    #     err = subprocess.call("echo 'alias load_cantera_mod=\"module purge && module load compiler/gcc/11.2.0\"' >> ~/.bashrc", shell=True)
    #     if not err :
    #         print("Alias load_cantera_mod added to ~/.bashrc")
    #     else :
    #         print("FAILED to add alias load_cantera_mod to ~/.bashrc, you can do it manually :\
    #               alias load_cantera_mod='module purge && module load compiler/gcc/11.2.0'")
    # else :
    #     print("Alias load_cantera_mod NOT added to ~/.bashrc")

if not error:
    print("*" * 80)
    print()
    print("To use this brandnew Cantera installation, you must update some of your environment variables.")
    print("To do so, you can add the following lines to your .bashrc (or equivalent):")
    text="""\
#cantera-avbp-3.1
function load_cantera
{7}
    source {4}

    {5}

    export PYTHONPATH={0}/lib/python{1}.{2}/site-packages:$PYTHONPATH 
    export PKG_CONFIG_PATH={0}/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH={0}/lib:$LD_LIBRARY_PATH
    export PATH={0}/bin:$PATH
    export PYTHON_CMD={3}

    #Only if you dont already have a custom lib folder:
    export CUSTOM_LIB={0}/mech_lib
    export LD_LIBRARY_PATH=$CUSTOM_LIB:$LD_LIBRARY_PATH

{6}
{8}
"""
    #get the output of which python to get the path to the python executable
    open_braket_string = "{"
    close_bracket_strinc = "}" 
    pythonpath = subprocess.check_output('which python3', shell=True).decode('utf-8').split('\n')[0]
    source_venv = subprocess.check_output('echo $VIRTUAL_ENV', shell=True).decode('utf-8').split('\n')[0] + "/bin/activate"
    if argument == 'NFS' and cluster == 'CALYPSO':
        module_gcc = "module load gcc/12.3.0"
    else:
        module_gcc = ""
    if argument == 'local':
        dyld_path = """\
    #Required for MacOS:
    export DYLD_LIBRARY_PATH=$CUSTOM_LIB
"""
    else:
        dyld_path = ""

    print()
    print(text.format(install_dir_path,str(python_version), str(sys.version_info[1]),pythonpath,source_venv,module_gcc,dyld_path,open_braket_string,close_bracket_strinc))
    print()
    print("DONT FORGET TO SOURCE YOUR .bashrc !")
    print()
    print("*" * 80)


README.txt           (docs directory)


Contents of the Directory
-------------------------------


   1 install_examples directory
     --------------------------------

      The install_examples directory contains various successful compilation environments that 
      Cantera has been tested out on. These files are "pre-preconfig files". In other words
      these are shell scripts that are executed before Cantera's main preconfig script is run.
      
      In order to use these scripts, you should copy them to the top directory of the Cantera
      distribution. Then, most of the scripts will still need some editing to specify path
      information and locations of the installation directory.
      Then, you can then execute the scripts to create the Makefiles in the distribution.
      A typical successful session of making the distribution, testing it, and then running
      the small test problem section would consist of the following commands:
 

         - Install a python math package (either numeric, numarray, or numpy)
         - Install sundials or choose vode within Cantera
         autoconf
         script_prepreconfig.sh
         make
         make install
         cd test_problems
         make 
         make test
       
      
      You can then go to the installation directories and test your Cantera installations
      further.


 

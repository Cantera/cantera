# -*- coding: utf-8 -*-

#################################################################################################
## Copyright (C) 2017 The University Of Dayton. All Rights Reserved.
##
## No part of this program may be photocopied, transferred, or otherwise reproduced in machine
## or human readable form without the prior written consent of The University Of Dayton.
#################################################################################################
#################################################################################################
## THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES,
## INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
## FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
## UNIVERSITY OF DAYTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
## INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
## OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
## NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
## EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  THE UNIVERSITY OF DAYTON
## HAS NO OBLIGATION TO SUPPORT THE SOFTWARE.
#################################################################################################

"""
Please cite this work as follows:
Briones, A.M., Olding, R., Sykes, J.P., Rankin, B.A., McDevitt, K., Heyne, J.S., 
"Combustion Modeling Software Development, Verification and Validation," Proceedings of the 2018 ASME Power & Energy Conference,
PowerEnergy2018-7433, Lake Buena Vista, FL. 
"""

import numpy as np
import sys
import getopt
import time
from collections import OrderedDict
from table_def import TableDef
from Utils import PrintFlush

def printhelp():
	# If the user requests help, display documentation.
    PrintFlush("""
    TableInversionFPV.py: Compute the inverted table from a Cantera flamelet table.

    Usage:
    TableInversionFPV.py [--input=<filename>]
                       [--output=<filename>]
                       [--help]
                       [--debug]

    Example:
    TableInversionFPV.py --input=flamelet.txt --output=newflamelet.txt

    An input file must be specified from the program to run. It should be of the FLUENT flamelet table format.
    All tables must have the same grid size for the program to work properly.	
	
    If the output file name is not given, one will be automatically generated.
    This automatic output file will be the same as the input file, with the name having "new" appended.

    """)
    sys.exit(1)

def tabletoarray(input, lastlocation):
    
	# Initialize the array to be constructed
    thearray=[]
    input.seek(lastlocation)
	
    while(1):
        try:
            # Read the next line and partition it into chunks divided by spaces
            lastposition = input.tell()
            templine = input.readline().rstrip().split(' ')
            for entry in templine:
                # Read each chunk as a float and add it to the end of the array
                temp = float(entry)
                thearray.append(temp)
        except Exception as e:
            # If a character cannot be parsed as a float, end the array
            input.seek(lastposition)
            break
    
	# Return the resulting array
    return thearray

def TableInvertFPV(options):	
    PrintFlush("Start flamelet Inversion of FPV type")

    timestart = time.clock()
        
    # read in table definition file
    table = TableDef()
    table.Read()
	
    """
    ********************* output control tab ****************************
    """
    YQOI = table.YQOI #user chooses species mass fraction to output
    RRQOI = table.RRQOI  # user chooses species net reaction rates (kg/m^3-s) to output
    """
    ***************************end of output control tab*******************************
    """

    if '--help' in options:
        printhelp()
	
    if '--input' in options:	
	    # If the user specifies an input file, open that file
        inputFile = options['--input']
        input = open(inputFile,'r')	
    else:
	    # Without an input specified, notify the user and close the program
        inputFile = None
        PrintFlush("You must specify an input file.")
        sys.exit(1)

    if '--output' in options:
	    # If the user specifies an output file, open that file
        outputFile = options['--output']
        file = open(outputFile,'w')
    else:
	    # Without an output file specified, create one automatically
        file = open("inverted_" + inputFile,'w')
	
	# Initialize variables that will be used during iteration through the input file
    tempstring = 'something'
    ZMean = []
    CMean = []	
    variable_dictionary = OrderedDict()
	
    while tempstring != '':	
		# Store the last read line in a convenient place	
        lastline = tempstring
        lastlocation = input.tell()
        tempstring = input.readline()
        try:
			# Attempt to read the next line and partition it into chunks divided by spaces
            templine = tempstring.rstrip().split(' ')
            for entry in templine:
                # Attempt to read each chunk as a float
                temp = float(entry)
            
			# If the reading is successful, store this chunk of text to an array
            temptable = tabletoarray(input, lastlocation)
            
			# If this array is also an independent variable, store it as such
            if lastline.startswith("MIXTURE"):
                if ZMean==[]:
                    # Read the first table of Cmean values, then ignore the rest
                    ZMean = temptable
                else:
                    ZMean = ZMean
            # Otherwise, attempt to add the array to the relevant variable's array
            else:
                try:
                    # If the variable is already in memory, append the new values
                    variable_dictionary[lastline.rstrip()].append(temptable)
                except Exception as e:
                    # If the variable is not already in memory, create an entry for it and start the new array
                    variable_dictionary[lastline.rstrip()] = [temptable]

        except Exception as e:
			# If the line cannot be parsed as a table, see if it can be parsed 
			# as a relevant flamelet parameter
            if templine[0].startswith("NUMOF"):
                numberofspecies = int(templine[-1])
            if templine[0].startswith("C") and not templine[0].startswith("Cp"):
                CMean.append(float(templine[-1]))
            if templine[0].startswith("PRESSURE"):
                pressure = float(templine[-1])

	# Debug line: make sure all keys are being captured from the input file
    # PrintFlush(list(variable_dictionary.keys()))

	# Reverse C and Z due to the normal operation of the diffusion domain
    ZMean.reverse()
    CMean.reverse()
    premix_chi = 0
	
    for mixture_fraction in ZMean:
	
        # Storing flamelet information
        # Header information for flamelets
        file.write('HEADER\n')
        file.write('PREMIX_CHI %.4f\n' % premix_chi)
        file.write('Z  %15.9e\n' % mixture_fraction)
        file.write('NUMOFSPECIES	%3d\n' % numberofspecies)
        file.write('GRIDPOINTS	%3d\n' % (len(CMean)))
        file.write('PRESSURE %6.2f\n' % pressure)
        file.write('BODY\n')
        file.write('REACTION_PROGRESS\n')

        # Write progress variable
        for j in range(0,len(CMean)):  
            progress_variable = CMean[j] #/ max(CMean)     # NON-NORMALIZED TABLE
            file.write('%15.9e '% progress_variable)
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')		
		
        # Write temperature
        file.write('\nTEMPERATURE\n')
        for j in range(0,len(CMean)):
            file.write('%15.9e '% variable_dictionary['TEMPERATURE'][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')			
       
    	# Write heat capacity (not possible in FLUENT)
        file.write('\nCp\n')
        for j in range(0,len(CMean)):
            file.write('%15.9e '% variable_dictionary['Cp'][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')	
			
	    # Write mean molecular weight (not possible in FLUENT)
        file.write('\nMMW\n')
        for j in range(0,len(CMean)):
            file.write('%15.9e '% variable_dictionary['MMW'][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')	

        # # Write density (not possible in FLUENT)
        file.write('\nrho\n')
        for j in range(0,len(CMean)):
            file.write('%15.9e '% variable_dictionary['rho'][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')
			
        # Write viscosity (not possible in FLUENT)
        file.write('\nmu\n')
        for j in range(0,len(CMean)):
            file.write('%15.9e '% variable_dictionary['mu'][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')
			
        # Write thermal conductivity (not possible in FLUENT)
        file.write('\nTC\n')
        for j in range(0,len(CMean)):
            file.write('%15.9e '% variable_dictionary['TC'][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')
			
    	# Write net production rates (not possible in FLUENT)
        for name in RRQOI:
            file.write('\nnet_rxn-'+name+'\n')
            for j in range(0,len(CMean)):
                file.write('%15.9e '% variable_dictionary['net_rxn-'+name][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
                if ((j+1)%5==0 and (j+1) != len(CMean)):
                    file.write('\n')				
    				
        # Write species mass fractions
        for name in YQOI:        # Sykes changed to include only select species in output
            species_name = 'massfraction-'+name
            file.write('\n'+species_name+'\n')
            for j in range(0,len(CMean)):
                file.write('%15.9e '% variable_dictionary[species_name][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
                if ((j+1)%5==0 and (j+1) != len(CMean)):
                    file.write('\n')
			
        # Writing reaction progress source (Needed for FLUENT FPI calculations)
        file.write('\nPREMIX_CDOT\n')
        file.write('%15.9e '% 0)
        for j in range(1,len(CMean)-1):
            file.write('%15.9e '% variable_dictionary['PREMIX_CDOT'][len(CMean)-j-1][len(ZMean) - ZMean.index(mixture_fraction) - 1])
            if ((j+1)%5==0 and (j+1) != len(CMean)):
                file.write('\n')
        file.write('%15.9e '% 0)

        # Generate a new premix_chi each time for the sake of FLUENT
        premix_chi += 0.0001
		
        file.write('\n')
        file.write('\n')    
    
    timeend = time.clock()
    PrintFlush('Table inversion total execution time (sec): ', (timeend-timestart))


if __name__ == "__main__":		
	
    # Declare the list of acceptable program options.
    longOptions = ['input=', 'output=', 'help', 'debug']

    # Parse the program arguments and place them in convenient format in "options"
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'dh', longOptions)
        options = dict()
        for o,a in optlist:
            options[o] = a
        if args:
            raise getopt.GetoptError('Unexpected command line option: ' + repr(' '.join(args)))
    except getopt.GetoptError as e:
        PrintFlush('Error parsing arguments:')
        PrintFlush(e)
        PrintFlush('Run "TableInversionFPV.py --help" to see options.')
        sys.exit(1)

    TableInvertFPV(options)
    sys.exit(0)
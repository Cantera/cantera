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
import os
#UDRI imports
from table_def import TableDef
from Utils import PrintFlush
from PDFScriptParallelCVarCppTrapz import FlameletsFileHandler
from MapProgressVar import ParallelRemeshing, RemeshTable, InputFileHandler
# UDRI C++ imports
from convolute import PyCreateDiffusionFlameletsFile

def printhelp():
	# If the user requests help, display documentation.
    PrintFlush("""
    TablePreMapFPV.py: Compute the inverted table from a Cantera flamelet table.

    Usage:
    TablePreMapFPV.py [--input=<filename>]
                       [--output=<filename>]
                       [--help]
                       [--debug]

    Example:
    TablePreMapFPV.py --input=flamelet.txt --output=newflamelet.txt

    An input file must be specified from the program to run. It should be of the FLUENT flamelet table format.
    All tables must have the same grid size for the program to work properly.	
	
    If the output file name is not given, one will be automatically generated.
    This automatic output file will be the same as the input file, with the name having "new" appended.

    """)
    sys.exit(1)
def WriteFlamelets(gridpoints, data, i, fast_file, keyList, pressure, numSpecies):

		# write file header
        fast_file.write_header(np.max(data[i*gridpoints:(i+1)*gridpoints,1]), numSpecies, gridpoints, pressure)

        for k, item in enumerate(keyList):
            out_array = np.ascontiguousarray(data[i*gridpoints:(i+1)*gridpoints,k], dtype = np.double)
            fast_file.write_data(out_array, item)

        #write out a few extra cr's
        fast_file.write_cr()
	
def TablePreMapFPV(options):	
    PrintFlush("Start flamelet Premapping for FPV table")
        
    # read in table definition file
    table = TableDef()
    table.Read()
	
    """
    ********************* output control tab ****************************
    """
    interpolation_type = table.interpolation_type
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
    else:
	    # Without an output file specified, create one automatically
        outputFile = 'premapped_' + inputFile
	
    # Read input file into a dictionary of species data
    inputFileHandler =  FlameletsFileHandler(inputFile);
    inputFileHandler.Read()
    variable_dictionary = inputFileHandler.variableDictionary
    ZMean = inputFileHandler.zMean
    CMean = inputFileHandler.cMean
    numSpecies = inputFileHandler.numSpecies
    pressure = inputFileHandler.pressure	

    #create list of keys with order rearranged
    keyList = list(variable_dictionary.keys())

    # create lengths
    zmeanLen = np.shape(variable_dictionary['MIXTURE_FRACTION'])[1]
    cmeanLen = np.shape(variable_dictionary['PROGRESS VARIABLE'])[0]
	
    header = 'FLAMELET PROGRESS VARIABLE\n\n' \
    +'[NUMBER_FMEAN]\n' \
    +str(zmeanLen) + '\n' \
    +'[NUMBER_CPROGRESS]\n' \
    +str(cmeanLen) + '\n' \
    +'[NUMBER_FVAR]\n' \
    +'1\n' \
    +'[NUMBER_CVAR]\n' \
    +'1\n' \
    +'[NUMBER_SPECIES]\n' \
    +str(len(table.YQOI)) + '\n' \
    +'[NUMBER_QOI]\n' \
    +str(len(keyList)-len(table.YQOI)) + '\n' \
    +'[NUMBER_VARIABLES]\n' \
    +str(len(keyList)+2) + '\n' \
    +'[DATA]\n' \
    +'fmean\t'+'cmean\t'+'fvar\t'+'cvar\t'

    for variable_name in keyList:
        if variable_name != 'MIXTURE_FRACTION' and variable_name != 'PROGRESS VARIABLE':
            header += variable_name + '\t'

    data = np.zeros((zmeanLen * cmeanLen, 4), dtype=np.double)
    for variable_name in keyList:
        if variable_name != 'MIXTURE_FRACTION' and variable_name != 'PROGRESS VARIABLE':
            data = np.hstack((data, np.reshape(variable_dictionary[variable_name], (zmeanLen * cmeanLen,1))))

    data[:,0:1]= np.reshape(variable_dictionary['MIXTURE_FRACTION'], (zmeanLen * cmeanLen,1))
    data[:,1:2]= np.reshape(variable_dictionary['PROGRESS VARIABLE'], (zmeanLen * cmeanLen,1))
		
    rdata = np.zeros(np.shape(data), dtype=np.double)

    for j in range(0, cmeanLen):
	    for i in range(0,zmeanLen):
			     rdata[j*zmeanLen + i,:] = data[j* zmeanLen + (zmeanLen - 1 - i), :] 

    for j in range(0, cmeanLen):
	    data[j* zmeanLen: (j+1)*zmeanLen,:] = rdata[(cmeanLen -1 - j)* zmeanLen: (cmeanLen - j)*zmeanLen, :] 
	
    #removing negative values of progress variable	
    for i in range(0, cmeanLen * zmeanLen):
        if data[i,1:2] < 1e-6 or i < zmeanLen:
            data[i,1:2] = 0.0

    input.close()
    tempfile = "temp.txt"
    np.savetxt(tempfile, data, fmt='%10.6e', header=header, comments='')
    inputFile = InputFileHandler(tempfile)
    inputFile.Read()

    B = RemeshTable(inputFile, interpolation_type, inputFile.data)
    os.remove(tempfile)

    #debug
    #np.savetxt(outputFile, B, fmt='%10.6e',header=header, comments='')

    B = np.delete(B, 2, axis=1) #remove fvar column
    B = np.delete(B, 2, axis=1) #remove cvar column
    B = np.flipud(B)	
    fast_file = PyCreateDiffusionFlameletsFile(outputFile)
    for i in range(zmeanLen):
        WriteFlamelets(zmeanLen, B, i, fast_file, keyList, pressure, numSpecies)
     
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
        PrintFlush('Run "TablePreMapFPV.py --help" to see options.')
        sys.exit(1)

    TablePreMapFPV(options)
    sys.exit(0)
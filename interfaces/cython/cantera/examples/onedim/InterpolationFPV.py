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

#python imports
import numpy as np
import sys
import getopt

from Utils import PrintFlush

def printhelp():
    # If the user requests help, display documentation.
    PrintFlush("""
Script for Alex.py: Find the maximum grid size in a FLUENT flamelet table.

Usage:
    Script for Alex.py [--input=<filename>]
                       [--output=<filename>]
                       [--help]
                       [--debug]

Example:
    Script for Alex.py --input=flamelet.txt --output=newflamelet.txt

An input file must be specified from the program to run. It should be of the FLUENT flamelet table format.    
    
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

# FPV interpolation function
def InterpolateFPV(options):

    PrintFlush("Start flamelet Interpolation of FPV type")

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
        output = open(outputFile,'w')
    else:
        # Without an output file specified, create one automatically
        output = open("New " + inputFile,'w')

    maximumGridSize = 0
    tempstring = 'something'

    # Cycle through all lines of the input file to find the maximum grid size
    while tempstring != '':
        # Read each line of the input file
        tempstring = input.readline()

        # When reaching a "gridpoints" line, isolate the numeric value
        if tempstring.startswith("GRID"):
            tempstring = ''.join(c for c in tempstring if c.isdigit())
            tempGridSize = int(tempstring)
        
            # Update the maximum grid size, if necessary
            if tempGridSize > maximumGridSize:
                maximumGridSize = tempGridSize
                input.readline()
                input.readline()
                input.readline()
                lastlocation = input.tell()
                maximumarray = tabletoarray(input, lastlocation)
                
    newX = maximumarray
    
    # Initialize variables to be used during file writing
    input.seek(0)
    tempstring = 'something'
    lastline = tempstring

    premix_chi = 0.0000
    
    # Cycle through all of the input file
    while tempstring != '':
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
            oldY = tabletoarray(input, lastlocation)

            if lastline.startswith("massfraction"):
                massfraction = True
            else:
                massfraction = False
            
            # If this array is also the independent variable, store it as such
            if lastline.startswith("MIXTURE"):
                oldX = oldY
                newY = newX
                #PrintFlush(oldX)
            else: 
                try:    
                    # Interpolate to get the new array    
                    #newY = interpolate(newX,oldX,oldY)
                    newY = np.flipud(np.interp(np.flipud(newX),np.flipud(oldX),np.flipud(oldY)))
                except Exception as e:
                    PrintFlush("Error in interpolation.",e)
            
            # Print the new array in the FLUENT format
            for i in range(0,len(newY)):
                if massfraction and newY[i] <= 0:
                    newY[i] = 0
                output.write('%15.9e '% newY[i])
                if ((i+1)%5==0 and (i+1) != len(newY)):
                    output.write('\n')
            output.write('\n')
            
        except Exception as e:
            # If the line cannot be parsed as a table, simply print it to the new file
            #PrintFlush(tempstring, e)
            if tempstring.startswith("GRID"):
                output.write('GRIDPOINTS    %3d\n' % (len(newX)))
            elif tempstring.startswith("PREMIX_CHI"):
                output.write('PREMIX_CHI ' + str(round(premix_chi,4)) + '\n')
                premix_chi += 0.0001
            else:
                output.write(tempstring)

    # write out conversion table
    WriteConversionTable(options['--output'])


def WriteConversionTable(file_name):
    file = open(file_name, "r")
    table = open("conversion_table.txt", "w")
    
    lambda_array = []
    templine = 'string'
    while templine != '':
        # Read each line of the output file
        templine = file.readline()

        # When reaching a "C" line, isolate the numeric value
        if templine.startswith("C"):
            templine = templine.rstrip().split('  ') 
            if (len(templine)>1):
                lambda_array.append(float(templine[1]))
    
    #back to the beginning of the file
    file.seek(0)
    
    mixture_fraction_array = []
    progress_variable_array = []
    with file as input_data:
        for line in input_data:
            if line.strip() == "MIXTURE_FRACTION":
                break
        for line in input_data:
            if line.strip() == "PROGRESS VARIABLE":
                break
            templine = line.rstrip().split(' ')
            for entry in templine:
                mixture_fraction_array.append(float(entry))
                
    c_size = len(mixture_fraction_array) * len(lambda_array)
    file = open(file_name, "r")
    with file as input_data:
        while len(progress_variable_array) < c_size:
            for line in input_data:
                if line.strip() == "PROGRESS VARIABLE":
                    break
            for line in input_data:
                if line.strip() == "TEMPERATURE":
                    break
                templine = line.rstrip().split(' ')
                for entry in templine:
                    progress_variable_array.append(float(entry))

    nGridPoints =  len(mixture_fraction_array)
    nFlamelets = len(lambda_array)
    print(nGridPoints, nFlamelets, len(progress_variable_array))
    for j in range(nFlamelets-1,-1,-1):
        for i in range(nGridPoints-1,-1,-1):
            table.write("%15.12f %15.12f %15.12f\n" % (abs(mixture_fraction_array[i]), \
			                             abs(progress_variable_array[nGridPoints* j+i]), \
			                              lambda_array[j]))

    table.close()
    
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
        PrintFlush('Run "Script for Alex.py --help" to see options.')
        sys.exit(1)

    InterpolateFPV(options)
    sys.exit(0)
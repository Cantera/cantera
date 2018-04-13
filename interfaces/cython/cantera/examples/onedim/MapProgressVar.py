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
## EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.? THE UNIVERSITY OF DAYTON
## HAS NO OBLIGATION TO SUPPORT THE SOFTWARE.
#################################################################################################

"""
Please cite this work as follows:
Briones, A.M., Olding, R., Sykes, J.P., Rankin, B.A., McDevitt, K., Heyne, J.S., 
"Combustion Modeling Software Development, Verification and Validation," Proceedings of the 2018 ASME Power & Energy Conference,
PowerEnergy2018-7433, Lake Buena Vista, FL. 
"""

import sys
import getopt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d,Axes3D
import tkinter
from tkinter import ttk
from multiprocessing import Process, Pipe
from collections import OrderedDict

# UDRI imports
from table_def import TableDef
from Utils import PrintFlush
from Graphics import PlotDlg

# Input FileHandler class - responsible for reading input file 
class InputFileHandler():

    def __init__(self, fileName):
        self._fileName = fileName
        self._numGridpoints = 0     # num fmean
        self._numFlamelets = 0      # num cmean
        self._numFvars = 0          # num fvar 
        self._numVars = 0           # num variables 
        self._numCvars = 1          # num of cvar
        self._numSpecies = 0        # num of species
        self._numQoi = 0            # num QOI
        self._columnHeader = "none" # column headings
        self._data = []             # table data
        self._input = open(fileName,'r')
        self._variableDict = OrderedDict()

    @property
    def numGridpoints(self):
        '''number of fmean values'''
        return self._numGridpoints

    @property
    def numFlamelets(self):
        '''number of cmean values'''
        return self._numFlamelets

    @property
    def numFvars(self):
        '''number of fvar values'''
        return self._numFvars

    @property
    def numVars(self):
        '''number of variables'''
        return self._numVars

    @property
    def numCvars(self):
        '''number of Cvar values'''
        return self._numCvars

    @property
    def numSpecies(self):
        '''number of species'''
        return self._numSpecies

    @property
    def numQoi(self):
        '''number of Qoi values'''
        return self._numQoi

    @property
    def columnHeader(self):
        '''column headings'''
        return self._columnHeader

    @property
    def data(self):
        '''table data'''
        return self._data


    # reads in a data set header and keeps track of number of lines read
    def _ReadHeader(self):
        count = 1
        done = False
        temp = self._input.readline()
        while done == False:
            name = temp.strip()
            if name == '[NUMBER_FMEAN]':
                temp = self._input.readline().strip()
                self._numGridpoints = int(temp)
                count = count + 1
            elif name == '[NUMBER_CPROGRESS]':
                temp = self._input.readline().strip()
                self._numFlamelets = int(temp)
                count = count + 1
            elif name == '[NUMBER_FVAR]':
                temp = self._input.readline().strip()
                self._numFvars = int(temp)
                count = count + 1
            elif name == '[NUMBER_CVAR]':
                temp = self._input.readline().strip()
                self._numCvars = int(temp)
                count = count + 1
            elif name == '[NUMBER_SPECIES]':
                temp = self._input.readline().strip()
                self._numSpecies = int(temp)
                count = count + 1
            elif name == '[NUMBER_QOI]':
                temp = self._input.readline().strip()
                self._numQoi = int(temp)
                count = count + 1
            elif name == '[NUMBER_VARIABLES]':
                temp = self._input.readline().strip()
                self._numVars = int(temp)
                count = count + 1
            elif name == '[DATA]':
                temp = self._input.readline().strip()
                self._columnHeader = temp
                done = True
                count = count + 1
            if done == False:
                temp = self._input.readline()
                count = count + 1
        return count

    # read in the input file header
    def Read(self):

        # read in the header and the data
        lineCount = self._ReadHeader()

        arrayLen = self._numCvars * self._numFvars * self._numFlamelets * self._numGridpoints
        arrayWidth = self._numVars
        self._data = np.zeros((arrayLen, arrayWidth), dtype=np.double)

        for i, line in enumerate(self._input):
            for j, temp in enumerate(line.split()):
                self._data[i,j] = float(temp)

        # close the file
        self._input.close()

# OutputFileHandler class - responsible for writing output file 
class OutputFileHandler():

    def __init__(self, fileName):
        self._fileName = fileName
        self._pressure = 0          # pressure
        self._numGridpoints = 0     # num fmean
        self._numFlamelets = 0      # num cmean
        self._numFvars = 0          # num fvar 
        self._numVars = 0           # num variables 
        self._numCvars = 1          # num of cvar
        self._numSpecies = 0        # num of species
        self._numQoi = 0            # num QOI
        self._columnHeader = "none" # column headings
        self._output = open(fileName,'w')

    @property
    def pressure(self):
        '''number of fmean values'''
        return self._pressure

    @pressure.setter
    def pressure(self, value):
        self._pressure = value

    @property
    def numGridpoints(self):
        '''number of fmean values'''
        return self._numGridpoints

    @numGridpoints.setter
    def numGridpoints(self, value):
        self._numGridpoints = value

    @property
    def numFlamelets(self):
        '''number of cmean values'''
        return self._numFlamelets

    @numFlamelets.setter
    def numFlamelets(self, value):
        self._numFlamelets = value

    @property
    def numFvars(self):
        '''number of fvar values'''
        return self._numFvars

    @numFvars.setter
    def numFvars(self, value):
        self._numFvars = value

    @property
    def numVars(self):
        '''number of variables'''
        return self._numVars

    @numVars.setter
    def numVars(self, value):
        self._numVars = value

    @property
    def numCvars(self):
        '''number of Cvar values'''
        return self._numCvars

    @numCvars.setter
    def numCvars(self, value):
        self._numCvars = value

    @property
    def numSpecies(self):
        '''number of species'''
        return self._numSpecies

    @numSpecies.setter
    def numSpecies(self, value):
        self._numSpecies = value

    @property
    def numQoi(self):
        '''number of Qoi values'''
        return self._numQoi

    @numQoi.setter
    def numQoi(self, value):
        self._numQoi = value

    @property
    def columnHeader(self):
        '''column headings'''
        return self._columnHeader

    @columnHeader.setter
    def columnHeader(self, value):
        self._columnHeader = value
      
    # write the header and body of file
    def Write(self, data):

        # Creating Header information for flamelets
        self._output.write("FLAMELET PROGRESS VARIABLE MAPPED\n\n")
        self._output.write("%s\n%.1f\n" % ("[PRESSURE]", self._pressure))
        self._output.write("%s\n%d\n" % ("[NUMBER_FMEAN]", self._numGridpoints))
        self._output.write("%s\n%d\n" % ("[NUMBER_CPROGRESS]", self._numFlamelets))
        self._output.write("%s\n%d\n" % ("[NUMBER_FVAR]", self._numFvars))
        self._output.write("%s\n%d\n" % ("[NUMBER_CVAR]", self._numCvars))
        self._output.write("%s\n%d\n" % ("[NUMBER_SPECIES]", self._numSpecies))
        self._output.write("%s\n%d\n" % ("[NUMBER_QOI]", self._numQoi))
        self._output.write("%s\n%d\n" % ("[NUMBER_VARIABLES]", self._numVars))
        self._output.write("%s\n%s\n" % ("[DATA]", self._columnHeader))

        arrayLen = self._numCvars * self._numFvars * self._numFlamelets * self._numGridpoints
        arrayWidth = self._numVars
        for i in range(0, arrayLen):
            for j in range(0, arrayWidth):
                self._output.write("%10.6e\t" % data[i, j])
            self._output.write("\n")

def ParallelTableMapping(xp_source, fp_source, x_target, f_target, interpolate, n_gridpoints,n_flamelets, outPipe):			
	
    for k in range(0,n_gridpoints):
        xp = []
        fp = []
        for i in range(0,n_flamelets):
            xp.append(xp_source[i*n_gridpoints+k])
            fp.append(fp_source[i*n_gridpoints+k])
        f =  interp1d(xp, fp, kind = interpolate, assume_sorted='False')
        for i in range(1,n_flamelets-1):
            x = x_target[i*n_gridpoints+k]
            f_target[i*n_gridpoints+k]=f(x)
    
    outPipe.send(f_target)
    outPipe.close()
    
			
def MapTable(fileHandler, interpolate):

    n_fvars = fileHandler.numFvars              # num fvars
    n_gridpoints = fileHandler.numGridpoints    # num fmean
    n_flamelets = fileHandler.numFlamelets      # num cmean
    n_cvars = fileHandler.numCvars              # num cvars

    # useful indexing constants
    n_gf = n_gridpoints * n_flamelets
    n_fgf = n_fvars * n_gridpoints * n_flamelets

    A = fileHandler.data                        # array of table data

    #the second column from A array is lambda
    A_Czst = np.copy(A[:,1])

    # read in conversion table
    B = np.loadtxt('conversion_table.txt')
			
    #the first column from B array is fmean
    B_fmean = np.copy(B[:,0])

    #the second column from B array is cmean
    B_cmean = np.copy(B[:,1])

    #the third column from B array is Czst
    B_Czst = np.copy(B[:,2])


    xp_source = B_Czst[:]
    fp_source = B_cmean[:]
    processDict = OrderedDict()
    for m in range(0, n_cvars):
        for n in range(0, n_fvars):
            x_target = A_Czst[m* n_fgf + n * n_gf:m* n_fgf +  (n+1) * n_gf]
            f_target = np.copy(x_target)
            inputPipe, outputPipe = Pipe(duplex=False)
            p = Process(target= ParallelTableMapping, args=(xp_source, fp_source, x_target, f_target, interpolate, n_gridpoints,n_flamelets, outputPipe))
            processDict[n+m*n_fvars] = [p, inputPipe]
            p.start()
		
    for m in range(0, n_cvars):
        for n in range(0, n_fvars):
             A[m* n_fgf + n * n_gf:m* n_fgf + (n+1) * n_gf,1] = processDict[n+m*n_fvars][1].recv()
	
    return A

def ParallelRemeshing(xp_source, fp_source, x_target, f_target, interpolate, n_cvars, n_fvars, n_flamelets, n_gridpoints, outPipe):
    # useful indexing constants
    n_gp2 = n_gridpoints * n_gridpoints
    n_fgp2 = n_fvars * n_gp2
    n_gf = n_gridpoints * n_flamelets
    n_fgf = n_fvars * n_gridpoints * n_flamelets

    for m in range(0, n_cvars):
        for n in range(0, n_fvars):
            for k in range(0,n_gridpoints):
                xp = []
                fp = []
                for i in range(0,n_flamelets):
                    xp.append(xp_source[m*n_fgf + n*n_gf + i*n_gridpoints+k])
                    fp.append(fp_source[m*n_fgf + n*n_gf + i*n_gridpoints+k])
                lower_bound = np.asscalar(fp[0])
                upper_bound = np.asscalar(fp[-1])
                bounds = (lower_bound, upper_bound)
                f = interp1d(xp,fp, kind = interpolate, bounds_error = False , fill_value = bounds, assume_sorted='False')
                for j in range(1,n_gridpoints):
                    x = x_target[m*n_fgp2 + n*n_gp2 + j*n_gridpoints+k]
                    f_target[m*n_fgp2 + n*n_gp2 + j*n_gridpoints+k] = f(x)

    outPipe.send(f_target)
    outPipe.close()	
						
def RemeshTable(fileHandler, interpolate, A):

    n_vars = fileHandler.numVars                # num variables
    n_gridpoints = fileHandler.numGridpoints    # fmean
    n_flamelets = fileHandler.numFlamelets      # num cmean
    n_fvars = fileHandler.numFvars              # num fvar
    n_cvars = fileHandler.numCvars              # num cvar

    # useful indexing constants
    n_gp2 = n_gridpoints * n_gridpoints
    n_fgp2 = n_fvars * n_gp2
    n_gf = n_gridpoints * n_flamelets
    n_fgf = n_fvars * n_gridpoints * n_flamelets	
	
    # allocate the B matrix
    n_mapped_matrix_size = n_fgp2 * n_cvars
    B = np.zeros((n_mapped_matrix_size , n_vars))

    #okay the table has been mapped in terms of cmean again, but the table is not order properly for the interpolation methods in Table.c
    #the mapped table needs to be re-interpolated on cmean again.

    #the maximum value of A[:,1] should be the maximum value of A_cmean
    max_A_cmean = np.max(A[:,1])

    # create a new range of cmean values based on the number of fmeans
    regular_grid_A_cmean = np.zeros((n_gridpoints,1))
    for i in range(0, n_gridpoints, 1):
        regular_grid_A_cmean[i] = i * max_A_cmean / (n_gridpoints - 1)

    fvar_values = np.zeros((n_fvars,1))	

    #the values of fvar are obtained
    for i in range(0, n_fvars, 1):
        fvar_values[i] = A[i * n_gf,2]

    cvar_values = np.zeros((n_cvars,1))	

    #the values of cvar are obtained
    for i in range(0, n_cvars, 1):
        cvar_values[i] = A[i * n_fgf,3]
   
    for i in range(0, n_gridpoints, 1):
        #copying values of fmean from A to B for fvar = 0 and cvar = 0
        B[i * n_gridpoints: (i + 1) * n_gridpoints, 0:1] = A[0:n_gridpoints,0:1]
        #copying values of cmean from regular_grid_A_cmean to B for fvar = 0 and cvar = 0
        B[i * n_gridpoints: (i + 1) * n_gridpoints, 1:2] = regular_grid_A_cmean[i:i+1,0:1]

    #copying values of fvar to B for cvar = 0		
    for i in range(0, n_fvars, 1):
        B[i * n_gp2: (i + 1) * n_gp2, 2:3] = fvar_values[i,0:1] 

    #copying values of fmean and cmean for all other fvar to B for cvar = 0
    for i in range(1, n_fvars, 1):
        B[i * n_gp2: (i+1) * n_gp2, 0:2] = B[0:n_gp2, 0:2]
		
    #copying values of cvar to B		
    for i in range(0, n_cvars, 1):
        B[i * n_fgp2: (i + 1) * n_fgp2, 3:4] = cvar_values[i, 0:1] 

    #copying values of fmean, cmean and fvar for all other cvar within B table
    for i in range(1, n_cvars, 1):
        B[i * n_fgp2: (i+1) * n_fgp2, 0:3] = B[0:n_fgp2, 0:3]

    # i in the old cmean index in A array
    # j is the new cmean index in B array		
    # k is the fmean index in either array (same dimensions)		
    # l is the variable index 
    # n is the fvar index
	# m is the cvar index

    # No need to interpolate the table for the lower branch at each cvar, fvar		
    for m in range(0, n_cvars):
        for n in range(0, n_fvars):
            B[m*n_fgp2 + n*n_gp2 : m*n_fgp2 + n*n_gp2 + n_gridpoints, 4:n_vars] = A[m*n_fgf + n*n_gf:m*n_fgf + n*n_gf + n_gridpoints,4:n_vars]
            B[m*n_fgp2 + (n+1)*n_gp2 - n_gridpoints: m*n_fgp2 + (n+1)*n_gp2,4:n_vars] = A[m*n_fgf + (n+1)*n_gf-n_gridpoints:m*n_fgf + (n+1)*n_gf,4:n_vars]

    processDict = OrderedDict()
    # multiple process parallel section
    for var in range(4,n_vars):
        xp_source = A[:,1]
        fp_source = A[:,var]
        x_target = B[:,1]
        f_target = B[:,var]
        inputPipe, outputPipe = Pipe(duplex=False)
        p = Process(target= ParallelRemeshing, args=(xp_source, fp_source, x_target, f_target, interpolate, n_cvars, n_fvars, n_flamelets, n_gridpoints, outputPipe))
        processDict[var] = [p, inputPipe]
        p.start()
		
    for var in range(4, n_vars):
         B[:,var] = processDict[var][1].recv()

    return B

# start of MapProgressVarStart function
def MapProgressVarStart(options):        

    PrintFlush("*** Map Progress Variable Start ***")

    # read in table definition file
    table = TableDef()
    if table.Read() != True:
        return
    
    # assign interpolation type
    interpolation_type = table.interpolation_type
    pressure = table.pressure

    PrintFlush("Successfully read table definition, interpolation = " + interpolation_type)

    if '--help' in options:
        #printhelp()
        placefiller = 1

    if '--input' in options:    
        # If the user specifies an input file, open that file
        inputFile = options['--input']
    else:
        # if not specified - use default
        inputFile = "UdriTable.fla"

    if '--output' in options:
        # If the user specifies an output file, open that file
        outputFile = options['--output']
    else:
        # Without an output file specified, create default automatically
        outputFile = "UdriTable.map"

    # read in input table fmean and lambda values
    PrintFlush("Reading input file - " + inputFile)
    inputFile = InputFileHandler(inputFile)
    inputFile.Read()

    # map the input file
    PrintFlush("Input file read complete - Mapping progress variable")
    A = MapTable(inputFile, interpolation_type)
    PrintFlush("Mapping progress variable complete - remeshing table")

    # remesh the table
    B = RemeshTable(inputFile, interpolation_type, A)
    PrintFlush("Remeshing progress variable complete - writing output file")

    # write the output file
    outputFile = OutputFileHandler(outputFile)
    outputFile.pressure = pressure
    outputFile.numGridpoints = inputFile.numGridpoints
    outputFile.numFlamelets = inputFile.numGridpoints
    outputFile.numFvars = inputFile.numFvars
    outputFile.numCvars = inputFile.numCvars
    outputFile.numQoi = inputFile.numQoi
    outputFile.numSpecies = inputFile.numSpecies
    outputFile.numVars = inputFile.numVars
    outputFile.columnHeader = inputFile.columnHeader
    outputFile.Write(B)
    PrintFlush("*** Progress Variable Recovered ***")
	
    # set plot mode to interactive
    plt.ion()
    init_plot = True
	
    # create lengths
    zmeanLen = inputFile.numGridpoints
    cmeanLen = inputFile.numGridpoints
    zvarLen = inputFile.numFvars
    cvarLen = inputFile.numCvars

    ZMean = np.zeros((zmeanLen,1))	
    #the values of fmean are obtained
    for i in range(0, zmeanLen, 1):
        ZMean[i] = B[i,0]	
		
    CMean = np.zeros((cmeanLen,1))	
    #the values of cmean are obtained
    for i in range(0, cmeanLen, 1):
        CMean[i] = B[i * zmeanLen,1]	

    ZVar = np.zeros((zvarLen,1))	
    #the values of fvar are obtained
    for i in range(0, zvarLen, 1):
        ZVar[i] = B[i * zmeanLen * cmeanLen,2]

    CVar = np.zeros((cvarLen,1))	
    #the values of cvar are obtained
    for i in range(0, cvarLen, 1):
        CVar[i] = B[i * zmeanLen * cmeanLen * zvarLen, 3]
	
    #create list of keys
    keyList = inputFile.columnHeader.split()[4:]
		
    # Initialize the dictionary containing variable matrix that will be used to visualize the output in 3D
    plotDict = OrderedDict()
    for varIndex, name in enumerate(keyList):
        plotDict[name] = np.zeros((cvarLen, zvarLen, cmeanLen, zmeanLen))
        for cvarIndex in range (cvarLen):           
            for zvarIndex in range(zvarLen):
                for cmeanIndex in range(cmeanLen):
                    for zmeanIndex in range(zmeanLen):
                         plotDict[name][cvarIndex, zvarIndex, cmeanIndex, zmeanIndex] = B[zmeanIndex + zmeanLen * cmeanIndex + zmeanLen * cmeanLen * zvarIndex + zmeanLen * cmeanLen * zvarLen * cvarIndex, 4+varIndex]
					
    # display plot dialog
    dlg = PlotDlg(CVar, ZVar, 0, 0, keyList, keyList[0])
    retrySolve = dlg.WaitResult()
    cvarIndex = dlg.cvar_index
    zvarIndex = dlg.zvar_index
    name = dlg.variable_name

    plotLabel = {}
    plotLabel['TEMPERATURE'] = 'Temperature, T [K]'
    plotLabel['mu'] = 'Dynamic Viscosity, $\mu$ [kg/m-s]'
    plotLabel['TC'] = 'Thermal Conductivity, k [W/m-K]'
    plotLabel['PREMIX_CDOT'] = 'Progress Variable Source, S$_C$ [kg/m^3-s]'
    plotLabel['MMW'] = 'Molecular Weight, MW [kg/kmol]'
    plotLabel['rho'] = r'$\mathrm{Density, \rho  [kg/m^3]}$'
    plotLabel['Cp'] = 'Specific Heat Capacity, C$_p$ [J/kg-K]'

    # plot variable profile
    while retrySolve == True:
        PrintFlush("plotting variable = ", name, " cvar = ", cvarIndex, " zvar = ", zvarIndex)
        if name in plotLabel:
            label = plotLabel[name]
        elif name.find('net_rxn-') != -1:
            label = name.replace('net_rxn-', '')
            label = label + ' Net Reaction Rate [kg/m$^3$-s]'
        elif name.find('massfraction-') != -1:
            label = name.replace('massfraction-', '')
            label = label + ' Mass Fraction'
        else:
            label = name
        fig = plt.figure(0)
        ax = Axes3D(fig)
        ZMeanAxis,CMeanAxis = np.meshgrid(ZMean,CMean)
        ax.plot_surface(ZMeanAxis,CMeanAxis,plotDict[name][cvarIndex, zvarIndex,:,:],rstride=1,cstride=1,cmap=cm.jet)
        plt.ylabel(r'Progress Variable, C')
        plt.xlabel("Mixture Fraction, Z")
        plt.xlim(0,1)
        plt.ylim(0,1)
        ax.set_zlabel(label)
        ax.view_init(15,-130)
        plt.title(str(table.fuel)+"-"+str(table.oxidizer)+", T$_{Fuel}$="\
        +str(round(table.Tfuel))+"K"+", T$_{Oxidizer}$="+str(round(table.Toxidizer))\
        +"K, P$_{Op}$="+str(round(table.pressure/101325,1))+"atm")
        txtstr = r'Scaled Z$_{var}$ = '+str(np.asscalar(ZVar[zvarIndex])) \
               +r'  Scaled C$_{var}$ = '+str(np.asscalar(CVar[cvarIndex]))
        ax.text2D(0.30, 0.87, txtstr, transform=ax.transAxes, fontsize=12)
        if init_plot:
            # only position initially, allows user to move afterwards
            mgr = plt.get_current_fig_manager()
            mgr.window.wm_geometry("+50+400")
            init_plot = False
       
        # display plot dialog
        dlg = PlotDlg(CVar, ZVar, cvarIndex, zvarIndex, keyList, name)
        retrySolve = dlg.WaitResult()        
        cvarIndex = dlg.cvar_index
        zvarIndex = dlg.zvar_index
        name = dlg.variable_name

# start of main
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
        PrintFlush('Run "MapProgressVar.py --help" to see options.')
        sys.exit(1)

    MapProgressVarStart(options)
    sys.exit()

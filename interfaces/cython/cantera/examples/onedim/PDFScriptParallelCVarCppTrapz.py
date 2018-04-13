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

# python imports
from scipy.stats import beta as Beta
from scipy.interpolate import RectBivariateSpline    # Potential use in Stanford flamelet models B, C, and D
from scipy.interpolate import UnivariateSpline
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys
import getopt
from mpl_toolkits.mplot3d import axes3d,Axes3D
import time
import tkinter
from tkinter import ttk
from multiprocessing import Process, Pipe
from collections import OrderedDict

# UDRI imports
from convolute import PyConvolute, PyCreateOutputFile
from table_def import TableDef
from Utils import PrintFlush
from Graphics import PlotDlg
   
# FlameletsFileHandler class - responsible for reading input file 
# and creating a dictionary of variable data 
# used by all processes
class FlameletsFileHandler():

    def __init__(self, fileName):
        self._numSpecies = 0
        self._zMean = []
        self._cMean = []
        self._pressure = 0
        self._numSpecies = 0 		
        self._gridPoints = 0
        self._variableDict = OrderedDict()
        self._input = open(fileName,'r')

    @property
    def variableDictionary(self):
        '''dictionary containing all species values'''
        return self._variableDict

    @property
    def zMean(self):
        '''zMean'''
        return self._zMean

    @property
    def cMean(self):
        '''cMean'''
        return self._cMean

    @property
    def pressure(self):
        '''pressure'''
        return self._pressure
		
    @property
    def numSpecies(self):
        '''numSpecies'''
        return self._numSpecies
		
    # reads in a data set header
    def _ReadHeader(self):
        body = False
        temp = self._input.readline()
        while body == False and temp != "":
            textArray = temp.strip().split()
            name = textArray[0]
            if name == 'NUMOFSPECIES':
                self._numSpecies = int(textArray[1])
            elif name == 'PRESSURE':
                self._pressure = float(textArray[1])
            elif name == 'Z':
                self._zMean.append(float(textArray[1]))
            elif name == 'GRIDPOINTS':
                self._gridPoints = int(textArray[1])
            elif name == 'BODY':
                body = True
            if body == False:
                temp = self._input.readline()
        return body

    # reads in all data associatted with a variable
    def _ReadVariable(self):
        varName = self._input.readline().strip()
        varArray = []
        if varName != "":
            count = 0
            while (count < self._gridPoints):
                templine = self._input.readline().strip().split(' ')
                count += len(templine)
                for entry in templine:
                    # convert text to float and add it to the end of the array
                    temp = float(entry)
                    varArray.append(temp)
        return varName, varArray

    # read in the input file and assign data to a dictionary of variables
    def Read(self):
        body = self._ReadHeader()
        while body == True:    
            varName, varArray = self._ReadVariable()
            while varName != "":
                if varName == 'REACTION_PROGRESS':
                    if self._cMean == []:
                        self._cMean = varArray
                else:
                    if varName in self._variableDict.keys():
                        self._variableDict[varName].append(varArray)
                    else:
                        self._variableDict[varName] = [varArray]
                varName, varArray = self._ReadVariable()
            body = self._ReadHeader()
        self._input.close()


# Determine betas for Z
def BetaPDF(realization, MEAN, GAMMA):
    alpha = MEAN*GAMMA
    beta = (1-MEAN)*GAMMA
    output = Beta.pdf(realization,alpha,beta)
    return output
    
def BetaCDF(realization, MEAN, GAMMA):
    alpha = MEAN*GAMMA
    beta = (1-MEAN)*GAMMA
    output = Beta.cdf(realization,alpha,beta)
    return output

def GetZVar(nZVar):
    ZVar = [0]      # Not used in laminar implementation
    if nZVar > 1:
        variances = np.logspace(-6.64386,0.0,nZVar-1,endpoint=True,base=2)
        for log in variances:
            ZVar.append(round(log,4))
    return ZVar

def GetCVar(nCVar):
    CVar = [0]      # Not used in Stanford implementation
    if nCVar > 1:
        variances = np.logspace(-6.64386,0.0,nCVar-1,endpoint=True,base=2)
        for log in variances:
            CVar.append(round(log,4))
    return CVar
    

# convolution - create as a process - one for each variable
def Convolution(name, filename, nZVar, nCVar, outPipe):

    PrintFlush("Convolution process active - initializing: " + name)

    # read input file into a dictionary of variable data
    inputFileHandler =  FlameletsFileHandler(filename)
    inputFileHandler.Read()
    variableDict = inputFileHandler.variableDictionary
    ZMean = inputFileHandler.zMean
    CMean = inputFileHandler.cMean
    zmeanLen = len(ZMean)
    cmeanLen = len(CMean)

    # start initialization timer
    betastart = time.clock()

    # Prepare convenience variables for the Farve-averaged calculations
    x = np.array(ZMean)
    y = np.array(CMean)
    xLen = len(x)
    yLen = len(y)
    max_y = max(y)

    # Initialize statistical properties of interest    
    ZVar = GetZVar(nZVar)
    CVar = GetCVar(nCVar)

    #convenient lengths of cvar and zvar
    zvarLen = len(ZVar)
    cvarLen = len(CVar)
    
    # calculate Gamma
    Gamma = np.zeros((len(ZVar)))
    for zvarIndex, zvar in enumerate(ZVar[1:]):
        if zvar != 1.0:
            Gamma[zvarIndex+1] = (1/zvar) - 1

    # calculate Gamma2
    Gamma2 = np.zeros((len(CVar)))
    for cvarIndex, cvar in enumerate(CVar[1:]):
        if cvar != 1.0:
            Gamma2[cvarIndex+1] = (1/cvar) - 1

    # Calculate values of the beta distribution
    Betas = np.zeros((zmeanLen, zmeanLen, zvarLen), dtype = np.double)
    for i in range(1,zvarLen):
        for j in range(1, zmeanLen-1):
            Betas[:, j, i] = BetaPDF(x,ZMean[j],Gamma[i])

    # Correcting the edges of the distribution from inside outwards
    for i in range(1,zvarLen):
        for j in range(1, zmeanLen-1):
            for k in range(xLen):
                if np.isinf(Betas[k,j,i]) or Betas[k,j,i] > 20.0:
                    Betas[k,j,i] = 20.0  
        
    Betas2 = np.zeros((cmeanLen, cmeanLen, cvarLen), dtype=np.double)
    for i in range(1, cvarLen):
        for j in range(1, cmeanLen-1):
            Betas2[:, j, i] = BetaPDF(y/max_y,CMean[j]/max_y,Gamma2[i])            
    
    # Correcting the edges of the distribution from inside outwards
    for i in range(1,cvarLen):
        for j in range(1, cmeanLen-1):
            for k in range(yLen):
                if np.isinf(Betas2[k,j,i]) or Betas2[k,j,i] > 20.0:
                    Betas2[k,j,i] = 20.0
    """
    if name == 'TEMPERATURE':
        fig = plt.figure()
        for i in range(1,zvarLen):
             for j in range(1, zmeanLen-1):
                plt.title("Zvar = %6.3f and Zmean = %12.8f" % (ZVar[i], ZMean[j]))
                plt.plot(x, Betas[:,j,i], figure = fig)
                plt.pause(1)
                plt.clf()
    """
    betaend = time.clock()
    PrintFlush('Initialization time (secs): ' + str((betaend-betastart)))

    # create output arrays
    zx = np.zeros(len(x))
    zy = np.zeros(len(y))
    zz = np.zeros((len(x), len(y)))

    # retrieve data from dictionary for this species
    dictEntry = np.array(variableDict[name])

    # create contiguous arrays for use by C++
    cBetas = np.ascontiguousarray(Betas, dtype=np.double)
    cBetas2 = np.ascontiguousarray(Betas2, dtype = np.double)
    cDictEntry = np.ascontiguousarray(dictEntry, dtype = np.double)
    cZz = np.ascontiguousarray(zz, dtype = np.double)
    cZx = np.ascontiguousarray(zx, dtype = np.double)
    cZy = np.ascontiguousarray(zy, dtype = np.double)

    # create C++ processing object and tell it what arrays to use
    conv = PyConvolute(name)
    conv.set_betas(cBetas)
    conv.set_betas2(cBetas2)
    conv.set_dict_entry(cDictEntry)
    conv.set_zz(cZz)
    conv.set_zx(cZx)
    conv.set_zy(cZy)

    # create array for zmean results to send to pipe
    pipe_results = np.zeros((zmeanLen), dtype=np.double)
    
    PrintFlush("Convolution process initialization complete for: " + name)

    maxValue =  np.max(dictEntry)    
    
    # process the input data for the given species 
    # and send array of zmean results to the file process        
    for cvarIndex, cvar in enumerate(CVar):
        for zvarIndex, zvar in enumerate(ZVar):
            for cmeanIndex, cmean in enumerate(CMean):
                for zmeanIndex, zmean in enumerate(ZMean):
                    if zvar == 0.0 and cvar == 0.0:
                        result = dictEntry[zmeanIndex, cmeanIndex]
                    elif zvar != 1.0 and cvar != 1.0:
                        conv.eval(cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex)
                    if zvar != 0.0 and zvar != 1.0 and cvar == 0.0:
                        result = trapz(cZx,x)
                    if zvar == 0.0 and cvar != 0.0 and cvar != 1.0:
                        result = trapz(cZy,y)/max_y
                    if zvar != 0.0 and zvar != 1.0 and cvar != 0.0 and cvar != 1.0:
                        result = trapz(trapz(Betas[:,zmeanIndex,zvarIndex] * dictEntry[zmeanIndex,cmeanIndex],x)*Betas2[:,cmeanIndex,cvarIndex], y/max_y)
                    if zvar == 1.0 or cvar == 1.0:
                        result = np.average(dictEntry[0,:]) * (1.0 - zmean) + np.average(dictEntry[-1,:]) * zmean                        
                    if np.isinf(result) or np.isnan(result):
                        result = 1.0e-30
                    if name == 'TEMPERATURE' or name == 'mu' or name == 'TC' or name == 'PREMIX_CDOT':
                        if result > maxValue:
                           result = maxValue

                    # save result in output array
                    pipe_results[zmeanIndex] = result
                
                
                leftBndry = dictEntry[0,cmeanIndex]    #air side
                rightBndry = dictEntry[-1,cmeanIndex]  #fuel side
                if name == 'rho':
                    leftBndry = 1 / leftBndry
                    rightBndry = 1 / rightBndry

                # Corrections to the edges of the tables are needed when convoluted
                if (zvar != 0.0 and zvar != 1.0) or (cvar != 0.0 and cvar != 1.0):
                    if name == 'TEMPERATURE' or name == 'mu' or name == 'TC' or name == 'Cp' or name == 'rho':
                        for zmeanIndex in range(len(ZMean)):
                            if pipe_results[zmeanIndex] < leftBndry:
                                pipe_results[zmeanIndex] = leftBndry * (1.0 - ZMean[zmeanIndex]) + rightBndry * ZMean[zmeanIndex]
                            else:
                                break
                        for zmeanIndex in range(len(ZMean)-1,-1,-1):
                            if pipe_results[zmeanIndex] < rightBndry:
                                pipe_results[zmeanIndex] = leftBndry * (1.0 - ZMean[zmeanIndex]) + rightBndry * ZMean[zmeanIndex]
                            else:    
                                break    
                    if name == 'rho':
                        for zmeanIndex in range(len(ZMean)):
                            pipe_results[zmeanIndex] = 1 / pipe_results[zmeanIndex]
                    
                    elif name == 'MMW':
                        if rightBndry > leftBndry:
                            for zmeanIndex in range(len(ZMean)):
                                if pipe_results[zmeanIndex] < leftBndry:
                                    pipe_results[zmeanIndex] = leftBndry * (1.0 - ZMean[zmeanIndex]) + rightBndry * ZMean[zmeanIndex]
                                else:
                                    break
                            for zmeanIndex in range(int(len(ZMean)/2)+1,len(ZMean)):                        
                                if pipe_results[zmeanIndex] < pipe_results[zmeanIndex-1] or pipe_results[zmeanIndex] > rightBndry:
                                    pipe_results[zmeanIndex] = leftBndry * (1.0 - ZMean[zmeanIndex]) + rightBndry * ZMean[zmeanIndex]
                        elif rightBndry < leftBndry: 
                            for zmeanIndex in range(int(len(ZMean)/2),-1,-1):
                                if pipe_results[zmeanIndex] < pipe_results[zmeanIndex+1] or pipe_results[zmeanIndex] > leftBndry:
                                    pipe_results[zmeanIndex] = leftBndry * (1.0 - ZMean[zmeanIndex]) + rightBndry * ZMean[zmeanIndex]                          
                            for zmeanIndex in range(len(ZMean)-1,-1,-1):                        
                                if pipe_results[zmeanIndex] < rightBndry:
                                    pipe_results[zmeanIndex] = leftBndry * (1.0 - ZMean[zmeanIndex]) + rightBndry * ZMean[zmeanIndex]
                                else:
                                    break      
                        elif rightBndry == leftBndry:
                            for zmeanIndex in range(len(ZMean)):
                                if pipe_results[zmeanIndex] < leftBndry:
                                    pipe_results[zmeanIndex] = leftBndry

                    pipe_results[0] = dictEntry[0, cmeanIndex]
                    pipe_results[-1] = dictEntry[-1, cmeanIndex]
                    if cmeanIndex == 0:
                        for zmeanIndex in range(len(ZMean)):                        
                            pipe_results[zmeanIndex] = dictEntry[zmeanIndex, cmeanIndex]

                # send output array to main file process via pipe
                outPipe.send(pipe_results)
                
        
    PrintFlush("Convolution process exiting " + name)

# output axis data (cmean, fmean, etc.) supporting a compressed file format
def OutputAxisData(output, axisName, axisData):
    output.write('[%s]\n' % (axisName))
    columnCount = 0
    for axisValue in axisData:
        output.write('%10.6e\t'% axisValue)
        columnCount += 1
        if columnCount >= 20:
            output.write('\n')
            columnCount = 0
    if columnCount != 0:
        # terminate partial line
        output.write('\n')

# help
def printhelp():
    # If the user requests help, display documentation.
    PrintFlush("""
    PDFScript.py: Compute the PDF table from a Cantera flamelet table.

    Usage:
    Script for Alex.py [--input=<filename>]
                       [--output=<filename>]
                       [--help]
                       [--debug]
                       [--compress]

    Example:
    PDFScript.py --input=flamelet.txt --output=newflamelet.txt

    An input file must be specified from the program to run. It should be of the FLUENT flamelet table format.
    All tables must have the same grid size for the program to work properly.    
    
    If the output file name is not given, one will be automatically generated.
    This automatic output file will be the same as the input file, with the name having "new" appended.

    """)
    sys.exit(1)

# start of PdfScript function
def PdfScriptStart(options):        
    
    # set plot mode to interactive
    plt.ion()
    init_plot = True

    #start the clock
    timestart = time.clock()

    # read in table definition file
    table = TableDef()
    if table.Read() != True:
        return
    
    """
    Obtained from Flamelet tab: Pass in the number of variances to use for Z (mixture fraction) and C (progress variable)
    """    
    nZVar = table.nZVar
    nCVar = table.nCVar
    """
    End variable passing; note that "0" variance number is the same as a laminar calculation
    """

    if '--help' in options:
        #printhelp()
        placefiller = 1

    if '--input' in options:    
        # If the user specifies an input file, open that file
        inputFile = options['--input']
    else:
        # if not specified - use default
        inputFile = "inverted_diffusion_flamelets.fla"

    if '--output' in options:
        # If the user specifies an output file, open that file
        outputFile = options['--output']
    else:
        # Without an output file specified, create default automatically
        outputFile = "CVar FPV "+str(table.fuel)+", "+str(round(table.pressure/101325,1))+"atm, "+str(round(table.Tfuel))+"K Inlet.fla"
    output = open(outputFile,'w')

    if '--compress' in options:
        compressed = True
        PrintFlush('Creating compressed output file')
    else:
        compressed = False
        PrintFlush('Creating original format output file')

    # Initialize statistical properties of interest    
    ZVar = GetZVar(nZVar)
    PrintFlush("SMix=",ZVar)
    PrintFlush()
    CVar = GetCVar(nCVar)            
    PrintFlush("SMixC=",CVar)
    PrintFlush()

    # Read input file into a dictionary of species data
    inputFileHandler =  FlameletsFileHandler(inputFile);
    inputFileHandler.Read()
    variable_dictionary = inputFileHandler.variableDictionary
    ZMean = inputFileHandler.zMean
    CMean = inputFileHandler.cMean

    # create lengths
    zmeanLen = len(ZMean)
    cmeanLen = len(CMean)
    cvarLen = len(CVar)
    zvarLen= len(ZVar)

    #create list of keys with order rearranged
    keyList = list(variable_dictionary.keys())
    keyList.remove('PREMIX_CDOT');
    index = keyList.index('TC');
    keyList.insert(index+1, 'PREMIX_CDOT');

    # Initialize the dictionary containing variable matrix that will be used to visualize the output in 3D
    plotDict = OrderedDict()
    for name in list(variable_dictionary.keys()):
        plotDict[name] = np.zeros((cvarLen, zvarLen, cmeanLen, zmeanLen))

    # Create and start a Convolution process for each variable name and add it to dictionary.
    # Each process is assigned a half-duplex pipe.
    processDict = OrderedDict()
    for name in keyList:
        inputPipe, outputPipe = Pipe(duplex=False)
        p = Process(target= Convolution, args=(name, inputFile, nZVar, nCVar, outputPipe))
        processDict[name] = [p, inputPipe]
        p.start()

    # Print the output file header
    if table.model == "FPI":
        output.write('FLAMELET PROLONGATION OF ILDM\n\n')
    else:
        output.write('FLAMELET PROGRESS VARIABLE\n\n')
    output.write('[NUMBER_FMEAN]\n')
    output.write(str(len(ZMean)) + '\n')
    output.write('[NUMBER_CPROGRESS]\n')
    output.write(str(len(CMean)) + '\n')
    output.write('[NUMBER_FVAR]\n')
    output.write(str(len(ZVar)) + '\n')
    output.write('[NUMBER_CVAR]\n')
    output.write(str(len(CVar)) + '\n')
    output.write('[NUMBER_SPECIES]\n')
    output.write(str(len(table.YQOI)) + '\n')
    output.write('[NUMBER_QOI]\n')
    output.write(str(len(keyList)-len(table.YQOI)) + '\n')
    output.write('[NUMBER_VARIABLES]\n')
    output.write(str(len(keyList)+4) + '\n')  # 4 statistical variables
    if (compressed == True):
        OutputAxisData(output, 'FMEAN', ZMean)
        OutputAxisData(output, 'CMEAN', CMean)
        OutputAxisData(output, 'FVAR', ZVar)
        OutputAxisData(output, 'CVAR', CVar)
    output.write('[DATA]\n')
    if compressed == False:
        output.write('fmean\tcmean\tfvar\tcvar\t')    
    for variable_name in keyList:
        output.write(variable_name + '\t')
    output.write('\n')
    output.close()

    # create contiguous numpy arrays needed by C++ CreateOutputFile
    outputArray = np.zeros((len(keyList), zmeanLen), dtype=np.double)
    cOutputArray = np.ascontiguousarray(outputArray, dtype = np.double)
    cZMean = np.ascontiguousarray(ZMean, dtype = np.double)
    c_output = PyCreateOutputFile(outputFile, cOutputArray, cZMean, compressed)

    # Cycle through all possible statistical situations
    for cvarIndex, cvar in enumerate(CVar):           
        for zvarIndex, zvar in enumerate(ZVar):
            cmeanstart = time.clock()
            for cmeanIndex, cmean in enumerate(CMean):

                # receive zmean arrays from each of the convolution processes
                for outIndex, variable_name in enumerate(keyList):
                    cOutputArray[outIndex] = processDict[variable_name][1].recv()
                    plotDict[variable_name][cvarIndex, zvarIndex, cmeanIndex] = cOutputArray[outIndex]

                #write array of zmean values to output file
                c_output.write(cvar, zvar, cmean)
            
            cmeanend = time.clock()
            PrintFlush("***** cvar = ", cvar, " zvar = ", zvar," time (secs): " + str((cmeanend-cmeanstart)))
          
    # Close the output file
    c_output.close()

    timeend = time.clock()
    PrintFlush('Total execution time (hr): ', (timeend-timestart)/3600)
       
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
    plotLabel['PREMIX_CDOT'] = 'Progress Variable Source, S$_C$ [kg/m^3/s]'
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
        ax.plot_surface(ZMeanAxis,CMeanAxis,plotDict[name][cvarIndex, zvarIndex,:,:],rstride=1,cstride=1,cmap=cm.jet,antialiased=False)
        plt.ylabel(r'Progress Parameter, $\lambda$')
        plt.xlabel("Mixture Fraction, Z")
        plt.xlim(0,1)
        plt.ylim(0,1)
        ax.set_zlabel(label)
        ax.view_init(15,-130)
        plt.title(str(table.fuel)+"-"+str(table.oxidizer)+", T$_{Fuel}$="\
        +str(round(table.Tfuel))+"K"+", T$_{Oxidizer}$="+str(round(table.Toxidizer))\
        +"K, P$_{Op}$="+str(round(table.pressure/101325,1))+"atm")
        txtstr = r'Scaled Z$_{var}$ = '+str(ZVar[zvarIndex]) \
               +r'  Scaled C$_{var}$ = '+str(CVar[cvarIndex])
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
    longOptions = ['input=', 'output=', 'help', 'debug', 'nocompress']

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
        PrintFlush('Run "PDFScript.py --help" to see options.')
        sys.exit(1)

    PdfScriptStart(options)
    sys.exit()
       
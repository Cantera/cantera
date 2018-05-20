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

import cantera as ct
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt
import warnings
import matplotlib.cbook
import tkinter
import getopt
import tempfile

from table_def import TableDef
from Utils import MessageDlg, PrintFlush


# Suppress matplotlib warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# UDRI C++ imports
from convolute import PyCreateDiffusionFlameletsFile

class RetryDlg():

    # constructor - display retry dialog box
    def __init__(self, OnePointControl, TwoPointControl, UserMaxTempPercentage, FuelSign, OxidSign, MINdT, MAXdT, dT, oxidizer_flux, cmeanSkip):
        self._dlg = tkinter.Tk()
        self._dlg.config(borderwidth=8)
        self._dlg.title("Flame Control")

        # local point control
        PointControl = False
        if OnePointControl == True:
            PointControl = 0
        if TwoPointControl == True:
            PointControl = 1

        # property definitions
        self._point_control = tkinter.IntVar(self._dlg, value=PointControl)
        self._max_temp = tkinter.DoubleVar(self._dlg, value=UserMaxTempPercentage*100.0)
        self._fuel_sign = tkinter.IntVar(self._dlg, value=FuelSign)
        self._oxidizer_sign = tkinter.IntVar(self._dlg, value=OxidSign)
        self._mindt = tkinter.DoubleVar(self._dlg, value=MINdT)
        self._maxdt = tkinter.DoubleVar(self._dlg, value=MAXdT)
        self._dt = tkinter.DoubleVar(self._dlg, value=dT)
        self._mdot_oxidizer = tkinter.DoubleVar(self._dlg, value=oxidizer_flux)
        self._cmean_skip = tkinter.IntVar(self._dlg, value=cmeanSkip)
        self._retry = False

        #one/two point control
        tkinter.Radiobutton(self._dlg, text="One-point Control", command = self.OnePointCallBack, variable = self._point_control, value=0).grid(row=2, column=1, sticky="w", pady=0, padx=2, columnspan=2)
        tkinter.Radiobutton(self._dlg, text="Two-point Control",  command = self.TwoPointCallBack,variable = self._point_control, value=1).grid(row=3, column=1, sticky="w", pady=0, padx=2, columnspan=2)

        #max temperature entry
        tkinter.Label(self._dlg, text="Max Temp %", padx=2).grid( row=4, column=1, sticky="w")
        tkinter.Entry(self._dlg, textvariable=self._max_temp, width=10).grid(row=4, column=2, sticky="w", pady=5, padx=2)

        # fuel increase/decrease
        tkinter.Radiobutton(self._dlg, text="Increase Fuel", variable = self._fuel_sign, value=1).grid(row=5, column=1, sticky="w", pady=0, padx=2, columnspan=2)
        tkinter.Radiobutton(self._dlg, text="Decrease Fuel", variable = self._fuel_sign, value=-1).grid(row=6, column=1, sticky="w", pady=0, padx=2, columnspan=2)

        # optional oxidizer flux
        self._mdot_label = tkinter.Label(self._dlg, text="Oxidizer Flux")
        self._mdot_label.grid( row=7, column=1, sticky="w")
        self._mdot_variable = tkinter.Entry(self._dlg, textvariable=self._mdot_oxidizer, width=10)
        self._mdot_variable.grid(row=7, column=2, sticky="w", pady=2, padx=2)

        #optional oxidizer increase/decrease
        self._increase_oxid = tkinter.Radiobutton(self._dlg, text="Increase Oxidizer", variable = self._oxidizer_sign, value=1)
        self._increase_oxid.grid(row=8, column=1, sticky="w", pady=0, padx=2, columnspan=2)
        self._decrease_oxid = tkinter.Radiobutton(self._dlg, text="Decrease Oxidizer", variable = self._oxidizer_sign, value=-1)
        self._decrease_oxid.grid(row=9, column=1, sticky="w", pady=0, padx=2, columnspan=2)
        tkinter.Label(self._dlg, text="Min dT").grid( row=10, column=1, sticky="w")

        # remaining mandatory entries
        tkinter.Entry(self._dlg, textvariable=self._mindt, width=10).grid(row=10, column=2, sticky="w", pady=2, padx=2)
        tkinter.Label(self._dlg, text="Max dT").grid( row=11, column=1, sticky="w")
        tkinter.Entry(self._dlg, textvariable=self._maxdt, width=10).grid(row=11, column=2, sticky="w", pady=2, padx=2)
        tkinter.Label(self._dlg, text="dT").grid( row=12, column=1, sticky="w")
        tkinter.Entry(self._dlg, textvariable=self._dt, width=10).grid(row=12, column=2, sticky="w", pady=2, padx=2)
        tkinter.Label(self._dlg, text="CMean Skip").grid( row=13, column=1, sticky="w")
        tkinter.Entry(self._dlg, textvariable=self._cmean_skip, width=10).grid(row=13, column=2, sticky="w", pady=2, padx=2)

        # retry and exit buttons
        tkinter.Button(self._dlg, text = "Retry", command = self.RetryCallBack, width = 6).grid( row=14, column=1, pady=8)
        tkinter.Button(self._dlg, text = "Exit", command = self.ExitCallBack, width = 6).grid( row=14, column=2, pady=8)

        # activate optional entries/radio buttons
        if OnePointControl == True:
            self.OnePointCallBack()
        if TwoPointControl == True:
            self.TwoPointCallBack()

    @property
    def one_point_control(self):
        '''point control selected'''
        return (self._point_control.get() == 0)

    @property
    def two_point_control(self):
        '''point control selected'''
        return (self._point_control.get() == 1)

    @property
    def max_temp(self):
        '''UserMaxTempPercentage'''
        return self._max_temp.get()/100.0

    @property
    def fuel_sign(self):
        '''FuelSign'''
        return self._fuel_sign.get()

    @property
    def oxidizer_sign(self):
        '''OxidSign'''
        return self._oxidizer_sign.get()

    @property
    def mindt(self):
        '''MINdT'''
        return self._mindt.get()

    @property
    def maxdt(self):
        '''MAXdT'''
        return self._maxdt.get()

    @property
    def dt(self):
        '''dT'''
        return self._dt.get()
    
    @property
    def mdot_oxidizer(self):
        '''mdot oxidizer'''
        return self._mdot_oxidizer.get()
    
    @property
    def cmean_skip(self):
        '''cmean skip'''
        return self._cmean_skip.get()
    
    # wait for user to hit button
    def WaitResult(self):
        self._dlg.mainloop()
        return self._retry

    # OnePoint radio button handler
    def OnePointCallBack(self):
        self._increase_oxid.grid_remove()
        self._decrease_oxid.grid_remove()
        self._mdot_label.grid()
        self._mdot_variable.grid()

    # TwoPoint radio button handler
    def TwoPointCallBack(self):
        self._increase_oxid.grid()
        self._decrease_oxid.grid()
        self._mdot_label.grid_remove()
        self._mdot_variable.grid_remove()

    # retry button handler
    def RetryCallBack(self):
        self._retry = True
        self._dlg.quit()
        self._dlg.destroy()

    # exit buttonhandler
    def ExitCallBack(self):
        self._dlg.quit()
        self._dlg.destroy()

# start of DiffusionFlame, set plot mode to interactive
def DiffusionFlame(options):

    # start of main, set plot mode to interactive
    init_plot = True
    PrintFlush("Start flamelet generation of FPV type")

    # Set up a temporary data directory for later reference
    tempdir = tempfile.TemporaryDirectory()
    data_directory = tempdir.name + "\\"
#    data_directory = r"C:\temp/"
#    if not os.path.exists(data_directory):
#        os.makedirs(data_directory)

    # read in table definition file
    table = TableDef()
    table.Read()
    table.Report()    

    #Create a summary file for later reference
    summary = open("Summary.txt", 'w')
    summary.write('     Amean       Amax       Tmax      TFfix      TOfix     FuelSign     OxSign       %loc    dT\n')
    fast_file = PyCreateDiffusionFlameletsFile('diffusion_flamelets.fla')

    """
    ******************** beginning of inputs ************************
    """

    """
    The following entry is obtained from the C# dialog models tab
    """
    transport_model = table.transport_model # this can be either mix or multi
    PrintFlush(transport_model, table.transport_model)
    UnityLewisNumber = False
    if transport_model == 'Unity':
        UnityLewisNumber = True
    transport_model = 'Mix'
    gas = ct.Solution(table.reaction_mechanism)
    flow = table.flow
    """
    end of models tab
    """

    """
    The values below are coming from the C# dialog boundaries tab.
    """
    pressure = table.pressure # pressure [Pa], information passed from C# dialog GUI
    Tfuel =  table.Tfuel # unburned gas temperature [K], information pass from C# dialog GUI
    Toxidizer = table.Toxidizer # unburned gas temperature [K], information pass from C# dialog GUI
    mole_fractions = table.mole_fractions # this information is passed from C# dialog GUI
    fuel = table.fuel  #the user chooses the fuel, this should be an array; this are the species which are greater than zero
                       # in the GUI
    fuel_fraction = table.fuel_fraction # this fraction is passed from the GUI
    oxidizer = table.oxidizer #the user chooses the oxidizer, this should be array; these are species > 0
    oxidizer_fraction = table.oxidizer_fraction #these are the mole fractions chosen by user
    init_fuel_flux = table.init_fuel_flux #kg/m^2-s, this is only a FPV input
    """
    end of boundary information
    """

    """
    The values below are coming from the C# dialog numerics tab.
    """
    initial_grid = np.linspace(table.start, table.end, table.num_points) # m, originally 20cm
    tol_ss = [table.ss_rel,table.ss_abs] # [rtol atol] for steady-state problem
    tol_ts = [table.tr_rel,table.tr_abs] # [rtol atol] for time stepping
    ss_age = table.ss_age  #Set the maximum number of times the Jacobian will be used before it must be re-evaluated by the steady solver.
    ts_age = table.ts_age  #Set the maximum number of times the Jacobian will be used before it must be re-evaluated by the transient solver. 
    stepsize = table.stepsize # initial time step in seconds
    nstep = table.nstep # sequence of integer numbers
    loglevel = table.loglevel # amount of diagnostic output (0 to 8)
    refine = table.refine  # 'True' to enable refinement, 'False' to disable
    temperature_limit_extinction =  table.temperature_limit_extinction # K  #lets add this parameter to the numerics tab
    ratio = table.ratio
    slope = table.slope
    curve = table.curve
    prune = table.prune
    """
    end of boundary information
    """

    """
    The values below are coming from the C# dialog flamelet tab.
    For now only one entry is being passed because we are developing the laminar flamelet
    """
    z_composition = ['H']
    composition = table.composition  # progress variable definition
    z_points = table.z_points #number of fmean. This is the same as 21 in inital_grid variable above
    cmean_skip = table.cmean_skip # number of cmean_skip for cmean entry in the GUI
    """
    end of flamelet information
    """

    """
    FPV setups requires an additional control tab named Flame Control
    This table does not need to be accessible for FPI
    """
    StrainEq = table.strain_rate
    IncludeFlameEquilibrium = False
    precalculations = table.precalculations
    UserStrainFactor = table.user_strain_factor
    UpperBranchHomotopicIter = table.zero_order_steps
    UserMaxTempPercentage = table.user_max_temp_percentage
    FuelSign = table.fuel_sign
    OxidSign = table.oxidizer_sign
    OnePointControl = table.one_point_control
    TwoPointControl = table.two_point_control
    dT = table.dT
    FailedSolutionsToTerminate = table.failed_solutions_to_terminate
    maxNumberOfFlames = table.max_num_flames
    N_opt = table.nopt
    MAXdT = table.max_dT
    MINdT = table.min_dT
    OxidizerMassFlux = table.oxidizer_mass_flux;

    """
    ending Flame Control tab
    """

    """
    ********************* output control tab ****************************
    """
    YQOI = table.YQOI #user chooses species mass fraction to output
    RRQOI = table.RRQOI  # user chooses species net reaction rates (kg/m^3-s) to output
    """
    ***************************end of output control tab*******************************
    """

    """
    ***********************end of inputs ****************************************
    """

    """******************************************************************************
    *****************The settings for the Cantera flamelet calculations are set here
    """

    def get_C():
        #start the clock
        timestart = time.clock()
        index_array = np.zeros((len(composition)), dtype=np.int)
        for i,name in enumerate(composition):
            index_array[i] = gas.species_index(name)
        C_array = np.zeros(len(f.grid), dtype=np.double)
        for j in range(0,len(f.grid)-1):
            C=0.0
            for i in index_array:
                C += f.Y[i,j]
            if j != 0 and C < max(C_array) and C > 0.5: 
                break
            C_array[j] = C
        timeend = time.clock()
        #PrintFlush('get_max_C execution time: ', (timeend-timestart))
        return C_array   

    def get_mixture_fraction():
        #start the clock
        timestart = time.clock()

        # create arrays needed
        mixture_fraction_array =  np.zeros((len(f.grid)), dtype=np.double)
        atom_end = np.zeros((len(z_composition), len(f.grid)), dtype=np.double)
        denominator = np.zeros((len(z_composition)), dtype=np.double)

        # convert dictionary to array for faster execution
        denominator_Z = 0.0
        for i, atom in enumerate(z_composition):
            atom_end[i] = f.elemental_mass_fraction(atom) 
            denominator_Z += f.elemental_mass_fraction(atom)[0] - f.elemental_mass_fraction(atom)[-1]   
                 
        # calculate mixture fraction
        for j in range(0,len(f.grid)):              
            numerator_Z = 0.0
            for i in range(len(z_composition)):
                numerator_Z += atom_end[i,j] - atom_end[i, -1]
            mixture_fraction = numerator_Z / (denominator_Z + 10**-30)
            mixture_fraction_array[j] = mixture_fraction
        timeend = time.clock()
        #PrintFlush('get_mixture_fraction execution time: ', (timeend-timestart))
        return mixture_fraction_array

    # write current solution to output file quickly using cython/c++
    def FastLaminarFlamelets(max_C, mixture_fraction_array):

        #start the clock
        timestart = time.clock()

        # write file header
        gridpoints = len(f.grid)
        fast_file.write_header(max_C, len(gas.Y), gridpoints, gas.P)

        # write mixture fraction
        mixture_array = np.ascontiguousarray(mixture_fraction_array, dtype = np.double)
        for j in range(0,gridpoints):  
            if j != 0 and mixture_array[j] == mixture_fraction_array[j-1]:
                mixture_array[j] -= mixture_array[j] * 10**-9 
        fast_file.write_data(mixture_array, 'MIXTURE_FRACTION')

        # write progress variable
        out_array = np.ascontiguousarray(C_array, dtype = np.double)
        fast_file.write_data(out_array, 'PROGRESS VARIABLE')

        # Write temperature
        out_array = np.ascontiguousarray(f.T, dtype = np.double)
        fast_file.write_data(out_array, 'TEMPERATURE')

        # Write heat capacity (not possible in FLUENT)
        out_array = np.ascontiguousarray(f.cp_mass, dtype = np.double)
        fast_file.write_data(out_array, 'Cp')

        # Write mean molecular weight (not possible in FLUENT)
        temp = np.zeros((gridpoints), dtype = np.double)
        out_array = np.ascontiguousarray(temp, dtype = np.double)
        for j in range(0,gridpoints):
            out_array[j] = np.dot(f.X[:,j],gas.molecular_weights[:])
        fast_file.write_data(out_array, 'MMW')

        # Write density (not possible in FLUENT)
        out_array = np.ascontiguousarray(f.density_mass, dtype = np.double)
        fast_file.write_data(out_array, 'rho')

        # Write viscosity (not possible in FLUENT)
        out_array = np.ascontiguousarray(f.viscosity, dtype = np.double)
        fast_file.write_data(out_array, 'mu')

        # Write thermal conductivity (not possible in FLUENT)
        out_array = np.ascontiguousarray(f.thermal_conductivity, dtype = np.double)
        fast_file.write_data(out_array, 'TC')
               
        # Write net production rates (if any)
        if len(RRQOI) > 0:
            index_array = np.zeros((len(RRQOI)), dtype=np.int)
            for i,name in enumerate(RRQOI):
                index_array[i] = gas.species_index(name)
            for n,i in enumerate(index_array):
                npr_array = np.ascontiguousarray(f.net_production_rates[i,:], dtype = np.double)
                fast_file.write_data_production(gridpoints, npr_array, gas.molecular_weights[i], 'net_rxn-'+ RRQOI[n])

        # Write species mass fractions (if any)
        if len(YQOI) > 0:
            index_array = np.zeros((len(YQOI)), dtype=np.int)
            for i,name in enumerate(YQOI):
                index_array[i] = gas.species_index(name)
            for n,i in enumerate(index_array):
                out_array = np.ascontiguousarray(f.Y[i,:], dtype = np.double)
                fast_file.write_data_min(len(out_array), out_array, "massfraction-"+YQOI[n], 1e-10)

        # Writing reaction progress source (PREMIX_CDOT)
        index_array = np.zeros((len(table.composition)), dtype=np.int)
        for i,name in enumerate(table.composition):
            index_array[i] = gas.species_index(name)    
        npr_array = np.ascontiguousarray(f.net_production_rates, dtype = np.double)
        ia_array = np.ascontiguousarray(index_array, dtype = np.int32)
#        if flow == "Laminar":
        mw_array = np.ascontiguousarray(gas.molecular_weights, dtype = np.double)
        fast_file.write_laminar_data_reaction(npr_array, mw_array, ia_array)
#        else:
#            t_array = np.ascontiguousarray(f.T, dtype = np.double)
#            fast_file.write_turbulent_data_reaction(npr_array, t_array, ia_array, pressure)      

        # write out a few extra cr's
        fast_file.write_cr()

        timeend = time.clock()
        PrintFlush('LaminarFlamelets execution time: ', (timeend-timestart))
    
    def iloc(percent):
        Ttarget = percent * np.max(f.T)
        for i in range(0,len(f.grid),1):
            if f.T[i] > Ttarget:
                z0 = i
                break
        for i in range(len(f.grid)-1,0,-1):
            if f.T[i] > Ttarget:
                z1 = i
                break
        return z0, z1

    def stepControl(ieval, N_opt, dT, MAXdT, MINdT):
        failed = False
        xi = float(N_opt) / (ieval + 1)
        if xi < 0.5:
            xi = 0.5
        elif xi > 2.0:
            xi = 2.0
        dT *= xi
        if dT > MAXdT:  
            dT = MAXdT
            failed = True
        elif dT < MINdT:
            dT = MINdT  
            failed = True
        return failed, dT
    
    def setFuelBoundary():
        reactants = ''
        for i, molecules in enumerate(fuel):
            reactants += str(fuel[i]) + ':' + str(fuel_fraction[i])
            if len(fuel) > 1 and i < len(fuel) - 1:
                reactants += ','
        if mole_fractions:
            f.fuel_inlet.X = reactants
            gas.TPX = Tfuel, target_pressure, reactants
        else:
            f.fuel_inlet.Y = reactants
            gas.TPY = Tfuel, target_pressure, reactants
    
    def setOxidBoundary():
        reactants = ''
        for i, molecules in enumerate(oxidizer):
            reactants += str(oxidizer[i]) + ':' + str(oxidizer_fraction[i])
            if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                reactants += ','
        if mole_fractions:
            f.oxidizer_inlet.X = reactants
            gas.TPX = Toxidizer, target_pressure, reactants
        else:
            f.oxidizer_inlet.Y = reactants
            gas.TPY = Toxidizer, target_pressure, reactants

    def molarFuelAirRatio():
        setFuelBoundary()
        setOxidBoundary()
        x = 0.0
        y = 0.0
        for j in range(0, gas.n_elements):
            if gas.element_names[j] == 'C':
                for i in range(gas.n_species):
                    x += f.fuel_inlet.X[i] * gas.n_atoms(gas.species_names[i], gas.element_names[j])
            if gas.element_names[j] == 'H':
                for i in range(gas.n_species):
                    y += f.fuel_inlet.X[i] * gas.n_atoms(gas.species_names[i], gas.element_names[j])
        stoic_num_oxid_moles = x + y / 4.0
        
        return stoic_num_oxid_moles

    def stoicFuelAirRatio():
        """
        Here we are assuming that fuel is only on the left and oxidizer on the right.
        This formulation is only for diffusion flamelets
        """
        setFuelBoundary()
        setOxidBoundary()

        numerator = np.dot(f.fuel_inlet.X, gas.molecular_weights)           
        denominator = np.dot (f.oxidizer_inlet.X, gas.molecular_weights)
        denominator *= molarFuelAirRatio() * np.sum(f.oxidizer_inlet.X) / f.oxidizer_inlet.X[f.flame.component_index('O2')-5]

        return numerator /  denominator

    def flameEquilbrium():
        Z = np.zeros((len(initial_grid),1))   
        Zstoic = stoicFuelAirRatio()
        Zfuel = 1.0
        Zoxid = 0.0
        
        modified_fuel_fraction = np.zeros(len(fuel_fraction), dtype= 'float')
        modified_fuel_fraction[:] = fuel_fraction
        
        for k in range(len(initial_grid)):
            Tinter = (Tfuel - Toxidizer) * (f.grid[k] - f.grid[-1]) / (f.grid[0] - f.grid[-1]) + Toxidizer
            Z[k] = Zfuel + (initial_grid[k] - table.start) * (Zstoic - Zfuel) / (0.5 * (table.end - table.start) - table.start)

            if Z[k] < Zstoic:
                Z[k] = Zoxid + (initial_grid[k] - table.end) * (Zstoic - Zoxid) / (0.5 * (table.end - table.start) - table.end)

            if Z[k] < 1.0:
                phi = np.asscalar(Z[k]) / (1.0 -  np.asscalar(Z[k])) / Zstoic
                   
                if phi > 0.0:
                    reactants = ''
                    shared_species = []
                    # if species is in both the fuel and oxidizer streams Cantera cannot take duplicates
                    for j, fuel_species in enumerate(fuel):
                        for i, oxidizer_species in enumerate(oxidizer):
                            if oxidizer_species == fuel_species:
                                shared_species.append(fuel_species)
                                temp = molarFuelAirRatio() * oxidizer_fraction[i] / oxidizer_fraction[0] / phi                             
                                modified_fuel_fraction[j] = fuel_fraction[j] + temp
                    for i, molecules in enumerate(fuel):
                        reactants += str(fuel[i]) + ':' + str(modified_fuel_fraction[i])
                        if len(fuel) > 1 and i < len(fuel) - 1:
                            reactants += ','                    
                        if len(shared_species) < len(oxidizer):
                            reactants += ','
                    for i, molecules in enumerate(oxidizer):
                        if molecules not in shared_species:
                            reactants += str(oxidizer[i]) + ':' \
                            + str( molarFuelAirRatio() * oxidizer_fraction[i] / \
                            oxidizer_fraction[0]/ phi)
                            if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                                reactants += ','

                    if mole_fractions:
                        gas.TPX = Tinter, pressure, reactants
                    else:
                        gas.TPY = Tinter, pressure, reactants
                    gas.equilibrate('HP')
                
                    T = gas.T
                    f.set_value(1, 'T', k, T)
                    for j in range(gas.n_species):
                        f.set_value(1, 5 + j, k, gas.Y[j])
                else:
                    f.set_value(1, 'T', k, Toxidizer)
                    for j in range(gas.n_species):
                        f.set_value(1, 5 + j, k, f.oxidizer_inlet.Y[j])
            else:
                f.set_value(1, 'T', k, Tfuel)
                for j in range(gas.n_species):
                    f.set_value(1, 5 + j, k, f.fuel_inlet.Y[j])
                    
    def limits():
        file =  open("equilibrium_table.txt", "w")
        Zstoic = stoicFuelAirRatio()
        Zfuel = 1.0
        Zoxid = 0.0
        Z = np.linspace(Zoxid, Zfuel, 101)
        
        modified_fuel_fraction = np.zeros(len(fuel_fraction), dtype= 'float')
        modified_fuel_fraction[:] = fuel_fraction
        for n in range(len(Z)):
            if Z[n] == 0.0:
            #maximum value of progress variable (C) as function of mixture fractiton(Z)
                file.write("%f %f\n" % (Zoxid, 0.0)) 

            elif Z[n] > 0.0 and Z[n] < 1.0:
                phi = np.asscalar(Z[n]) / (1.0 -  np.asscalar(Z[n])) / Zstoic
                
                reactants = ''
                shared_species = []
                # if species is in both the fuel and oxidizer streams Cantera cannot take duplicates
                for j, fuel_species in enumerate(fuel):
                    for i, oxidizer_species in enumerate(oxidizer):
                        if oxidizer_species == fuel_species:
                            shared_species.append(fuel_species)
                            temp = molarFuelAirRatio() * oxidizer_fraction[i] / oxidizer_fraction[0] / phi                             
                            modified_fuel_fraction[j] = fuel_fraction[j] + temp
                for i, molecules in enumerate(fuel):
                    reactants += str(fuel[i]) + ':' + str(modified_fuel_fraction[i])
                    if len(fuel) > 1 and i < len(fuel) - 1:
                        reactants += ','                    
                    if len(shared_species) < len(oxidizer):
                        reactants += ','
                for i, molecules in enumerate(oxidizer):
                    if molecules not in shared_species:
                        reactants += str(oxidizer[i]) + ':' \
                        + str( molarFuelAirRatio() * oxidizer_fraction[i] / \
                        oxidizer_fraction[0]/ phi)
                        if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                            reactants += ','    
                if mole_fractions:
                    gas.TPX = 0.5*(Tfuel+Toxidizer), pressure, reactants
                else:
                    gas.TPY = 0.5*(Tfuel+Toxidizer), pressure, reactants
                gas.equilibrate('HP')
                
                C_limit = 0.0
                for species in composition:
                    C_limit += gas.Y[gas.species_index(species)]
                file.write("%f %f\n" % (Z[n], C_limit)) 

            elif Z[n] == 1.0:
            #mixture fraction(Z) vs. maximum value of progress variable (C)
                file.write("%f %f\n" % (Zfuel, 0.0)) 

        file.close()  
                       
    # Exponents for the initial solution variation with changes in strain rate
    # Taken from Fiala and Sattelmayer (2014)
    exp_d_a = - 1. / 2.
    exp_u_a = 1. / 2.
    exp_V_a = 1.
    exp_lam_a = 2.
    exp_mdot_a = 1. / 2.

    #Flame object
    f = ct.CounterflowDiffusionFlame(gas, initial_grid)
    f.set_max_grid_points(1, 250)
    f.set_max_jac_age(ss_age,ts_age)  #for some reason this was commented
    f.set_time_step(stepsize,nstep)  # this was not here 
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
    f.set_grid_min(1e-5) # we can keep this input to python script only
    f.transport_model = transport_model

    f.fuel_inlet.T = Tfuel
    f.oxidizer_inlet.T = Toxidizer
    f.fuel_inlet.spread_rate = 0.0
    f.oxidizer_inlet.spread_rate = 0.0
    target_pressure = pressure
    f.fuel_inlet.mdot = init_fuel_flux
    refine_grid = True # There is no need to add this one again
    if OnePointControl == True:
        oxidizer_flux = OxidizerMassFlux
    else:
        oxidizer_flux = f.oxidizer_inlet.mdot # oxidizer flux initialization

    #############################
    file_name = 'restart.xml'  
    restart = False
    retrySolve = True

    #############################################
    middleBranch = False
    lowerBranch = False
    turning = False
    flameControl = False
    # loop trying to solve for current temperatures
    # allows retry if a solution fails
    limits()
    while retrySolve == True:
        Continuation = True
        stepSizeControl = False # Active/Deactive step size control
        if restart == False:
            #first time in solve loop
            T_max = []
            a_max = []
            m_dots = []
       
            if IncludeFlameEquilibrium == True:
                flameEquilbrium()
                C_array = get_C()
                max_C = np.max(C_array)
                mixture_fraction_array = get_mixture_fraction()
                FastLaminarFlamelets(max_C, mixture_fraction_array)
                PrintFlush(-1, max(f.T), 0.000, "upper branch", "False")
                T_max.append(max(f.T))
                a_max.append(0.01)
                m_dots.append(0.0)
                if len(T_max) == 1:
                    global_max_C = max_C
                    iFlamesInUpperBranch = 1
                    iFlamesInMiddleBranch = 0
                    iFlamesInLowerBranch = 0
                    n = 1
                fig = plt.figure(0)
                fig.suptitle("Fuel-Oxidizer S-Curve")
                title_text = '{}-{} T{}={:6.1f}K, T{}={:6.1f}K, P{}={:6.0f}Pa'.format(fuel, oxidizer, "$_{Fuel}$", Tfuel, "$_{Oxidizer}$", Toxidizer, "$_{Op}$",pressure)
                plt.title(title_text, fontsize=10)
                plt.ylabel("Maximum Temperature, T$_{max}$ [K]")
                plt.xlabel("Strain Rate [1/s]")
                txtstr =  "$\Lambda$ = %6.4f\nT$_{F,fix}$=%6.1f\nT$_{Ox,fix}$=%6.1f" % (max_C, 0.0, 0.0)
                props = dict(boxstyle='round', facecolor='wheat', alpha=1.0)
                txt = plt.text(0.15, 0.60, txtstr, fontsize = 14, transform = fig.transFigure, verticalalignment = 'top', bbox = props)
                plt.semilogx(a_max, T_max, 'b', figure = fig)
                txt.set_visible(True)
                plt.pause(0.05)
                txt.set_visible(False)
                if init_plot:
                # only position initially, allows user to move afterwards
                     mgr = plt.get_current_fig_manager()
                     mgr.window.wm_geometry("+50+400")
                fig1 = plt.figure(1)
                fig1.suptitle("Counterflow Flames")
                plt.title(title_text, fontsize=10)
                plt.ylabel("Temperature, T [K]")
                plt.xlabel("Mixture Fraction, Z")
                plt.plot(mixture_fraction_array, f.T, figure = fig1)
                if init_plot:
                # only position initially, allows user to move afterwards
                     mgr = plt.get_current_fig_manager()
                     mgr.window.wm_geometry("+750+400")
                     init_plot = False
        
            for i in range(precalculations):
                PrintFlush(" Precalculations")
                target_pressure = pressure * (i + 1) / precalculations
                f.P = target_pressure
        
                reactants = ''
                for i, molecules in enumerate(fuel):
                    reactants += str(fuel[i]) + ':' + str(fuel_fraction[i])
                    if len(fuel) > 1 and i < len(fuel) - 1:
                        reactants += ','
                if mole_fractions:
                    f.fuel_inlet.X = reactants
                    gas.TPX = Tfuel, target_pressure, reactants
                else:
                    f.fuel_inlet.Y = reactants
                    gas.TPY = Tfuel, target_pressure, reactants
                fuel_density = gas.density_mass
        
                reactants = ''
                for i, molecules in enumerate(oxidizer):
                    reactants += str(oxidizer[i]) + ':' + str(oxidizer_fraction[i])
                    if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                        reactants += ','
                if mole_fractions:
                    f.oxidizer_inlet.X = reactants
                    gas.TPX = Toxidizer, target_pressure, reactants
                else:
                    f.oxidizer_inlet.Y = reactants
                    gas.TPY = Toxidizer, target_pressure, reactants
                oxidizer_density = gas.density_mass
        
                if i > 0:
                    f.fuel_inlet.mdot *= 1.0 # pow((i + 1)/i , 1.25)
                f.oxidizer_inlet.mdot = f.fuel_inlet.mdot * np.sqrt(oxidizer_density/fuel_density)
       
                ###########################################################################
                f.set_flame_control(1, StrainEq, UnityLewisNumber, False, False, -2000, -20, -2000, -20, False)
                if i == 0:
                    f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=True)
                else: 
                    f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=False)
                PrintFlush("  Pressure      Tmax Fuel_mdot")
                PrintFlush("%10.1f %10.1f %10.2f" % (f.P, max(f.T), f.fuel_inlet.mdot))
            
            f.save(data_directory + file_name, name='solution',
            description='Cantera version ' + ct.__version__ +
            ', reaction mechanism ' + table.reaction_mechanism, loglevel=loglevel)  
            
            mixture_fraction_array = get_mixture_fraction()	
            title_text = '{}-{} T{}={:6.1f}K, T{}={:6.1f}K, P{}={:6.0f}Pa'.format(fuel, oxidizer, "$_{Fuel}$", Tfuel, "$_{Oxidizer}$", Toxidizer, "$_{Op}$",pressure)
            fig1 = plt.figure(1)
            fig1.suptitle("Counterflow Flames")
            plt.title(title_text, fontsize=10)
            plt.ylabel("Temperature, T [K]")
            plt.xlabel("Mixture Fraction, Z")
            plt.plot(mixture_fraction_array, f.T, figure = fig1)
    

        else:
            # restarting solve loop after a failure
            flameControl = True
            middleBranch = True
            sign = -1
            f.restore(data_directory + file_name, name='solution', loglevel=loglevel)

        # Write the limits to define accessible regions in the table
        # Lists to save data
        T_max.append(max(f.T))
        a_max.append(f.strain_rate('max'))
        m_dots.append(f.fuel_inlet.mdot)

        C_array = get_C()
        max_C = np.max(C_array)
        if len(T_max) == 1:
            global_max_C = max_C
            iFlamesInUpperBranch = 1
            iFlamesInMiddleBranch = 0
            iFlamesInLowerBranch = 0
            n = 1

        mixture_fraction_array = get_mixture_fraction()
        FastLaminarFlamelets(max_C, mixture_fraction_array)

        iFailedSolutions = 0

        if iFlamesInUpperBranch == 1:
            fig = plt.figure(0)
            fig.suptitle("Fuel-Oxidizer S-Curve")
            title_text = '{}-{} T{}={:6.1f}K, T{}={:6.1f}K, P{}={:6.0f}Pa'.format(fuel, oxidizer, "$_{Fuel}$", Tfuel, "$_{Oxidizer}$", Toxidizer, "$_{Op}$",pressure)
            plt.title(title_text, fontsize=10)
            plt.ylabel("Maximum Temperature, T$_{max}$ [K]")
            plt.xlabel("Maximum Strain Rate, a$_{max}$ [1/s]")
            txtstr =  "$\Lambda$ = %6.4f\nT$_{F,fix}$=%6.1f\nT$_{Ox,fix}$=%6.1f" % (max_C, 0.0, 0.0)
            props = dict(boxstyle='round', facecolor='wheat', alpha=1.0)
            txt = plt.text(0.15, 0.60, txtstr, fontsize = 14, transform = fig.transFigure, verticalalignment = 'top', bbox = props)
            plt.semilogx(a_max, T_max, 'b', figure = fig)
            txt.set_visible(True)
            plt.pause(0.05)
            txt.set_visible(False)
            if init_plot:
            # only position initially, allows user to move afterwards
                 mgr = plt.get_current_fig_manager()
                 mgr.window.wm_geometry("+50+400")
            fig1 = plt.figure(1)
            fig1.suptitle("Counterflow Flames")
            plt.title(title_text, fontsize=10)
            plt.ylabel("Temperature, T [K]")
            plt.xlabel("Mixture Fraction, Z")
            plt.plot(mixture_fraction_array, f.T, figure = fig1)
            if init_plot:
            # only position initially, allows user to move afterwards
                 mgr = plt.get_current_fig_manager()
                 mgr.window.wm_geometry("+750+400")
                 init_plot = False
    
        iflag = 0
        f.clear_stats()
        while Continuation:
            stepFailed = False
            Tfuel_j, Toxid_j = iloc(UserMaxTempPercentage)
            Tfuel_fixed = f.T[Tfuel_j] + FuelSign * dT
            Toxid_fixed = f.T[Toxid_j] + OxidSign * dT
        
            if restart == False and iFlamesInUpperBranch > UpperBranchHomotopicIter:
                flameControl = True
        
            if n > 50 and turning == False:
                if np.sign((a_max[-1] - a_max[-2]) * (a_max[-2] - a_max[-3])) < 0:
                    middleBranch = True
   
            if middleBranch == True and iFlamesInMiddleBranch > 2000:
                strain_factor = a_max[-1] / a_max[-2]
                f.flame.grid *= strain_factor ** exp_d_a   
                                    
            if flameControl == False and middleBranch == False:
                strain_factor = (n + UserStrainFactor ) / n
                # Create an initial guess based on the previous solution
                # Update grid
                f.flame.grid *= strain_factor ** exp_d_a
                normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
                # Update mass fluxes
                f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
                f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
                # Update velocities
                f.set_profile('u', normalized_grid, f.u * strain_factor ** exp_u_a)
                f.set_profile('V', normalized_grid, f.V * strain_factor ** exp_V_a)
                # Update pressure curvature
                f.set_profile('lambda', normalized_grid, f.L * strain_factor ** exp_lam_a)

            if flameControl == True:
                oxidizer_flux = f.oxidizer_inlet.mdot
                if OnePointControl == True:
                    f.oxidizer_inlet.mdot = oxidizer_flux
                f.set_flame_control(1, StrainEq, UnityLewisNumber, OnePointControl, TwoPointControl, \
                                    Tfuel_fixed, Tfuel_j, Toxid_fixed, Toxid_j, True)

            try:
                success = True
                f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=False)  
        
            except Exception:
                txt.set_visible(True)
                success = False

            if success == True:
                txt.set_visible(False)
                if f.extinct() == False:
                    if n % 1 == 0:
                        T_max.append(max(f.T))
                        a_max.append(f.strain_rate('max'))
                        m_dots.append(f.fuel_inlet.mdot)
                    
                        summary.write("%10.1f %10.1f %10.1f %10.1f %10.1f %10d %10d %10.1f %10.1f\n" % (f.strain_rate('mean'), \
                        f.strain_rate('max'), max(f.T), Tfuel_fixed, Toxid_fixed, FuelSign, OxidSign, UserMaxTempPercentage, dT))
                
                    if n % 1 == 0 and T_max[-1] < T_max[-2]:
                        f.save(data_directory + file_name, name='solution',
                        description='Cantera version ' + ct.__version__ +
                        ', reaction mechanism ' + table.reaction_mechanism, loglevel=0)

                    C_array = get_C()
                    max_C = np.max(C_array)
                    mixture_fraction_array = get_mixture_fraction()
    
                    if n % cmean_skip == 0:
                        FastLaminarFlamelets(max_C, mixture_fraction_array)
                    
                    ieval = f.eval_count_stats[-1]
                    if middleBranch == False and lowerBranch == False:
                        branch = 'upper branch'
                    elif middleBranch == True and lowerBranch == False:
                        branch = 'middle branch'
                    elif middleBranch == False and lowerBranch == True:
                        branch = 'lower branch'
                    if n ==1 or n % 10 == 0:
                        PrintFlush("     #n      Tmax       Amax         dT     branch   FlameControl?   # gridpoints")
                    PrintFlush("%6d %10.1f %10.1f %10.2f %10s %10s %12d" \
                            % (n, max(f.T), f.strain_rate('max'), dT, branch.split(' ')[0], str(flameControl),  f.flame.n_points))
                    if middleBranch == False:
                        iFlamesInUpperBranch += 1
                    if middleBranch == True:
                        iFlamesInMiddleBranch += 1     
                        if T_max[-1] > T_max[-2]:
                            PrintFlush("***********Maximum temperature higher than that of previous solution******")
                            txt.set_visible(True)
                            T_max.pop(-1)
                            a_max.pop(-1)
                            m_dots.pop(-1)
                            Continuation = False
                        else:
                             summary.write("%10.1f %10.1f %10.1f %10.1f %10.1f %10d %10d %10.1f %10.1f\n" % (f.strain_rate('mean'), \
                             f.strain_rate('max'), max(f.T), Tfuel_fixed, Toxid_fixed, FuelSign, OxidSign, UserMaxTempPercentage, dT))
                    if iFailedSolutions > FailedSolutionsToTerminate:
                        Continuation = False
                else:          
                    if turning == False:
                        iflag = 1
                        f.restore(data_directory + file_name, name='solution', loglevel=loglevel)
                        flameControl = True
            else:
                txt.set_visible(True)
                iflag = 1
                iFailedSolutions += 1
                stepSizeControl = True     
                flameControl = True
                stepFailed, dT = stepControl(ieval, N_opt, dT, MAXdT, MINdT)
                PrintFlush('Calculation failed; the new value of dT is %f' % dT)
                f.restore(data_directory + file_name, name='solution', loglevel=loglevel)
                if stepFailed == False:
                    continue
                else:
                    Continuation = False

            if Continuation == True:
                # success - create or update plot, update number of flames, perform step control
                if success == True and T_max[-1] < T_max[-2]:
                    fig = plt.figure(0)
                    title_text = '{}-{} T{}={:6.1f}K, T{}={:6.1f}K, P{}={:6.0f}Pa'.format(fuel, oxidizer, "$_{Fuel}$", Tfuel, "$_{Oxidizer}$", Toxidizer, "$_{Op}$",pressure)
                    txtstr =  "$\Lambda$ = %6.4f\nT$_{F,fix}$=%6.1f\nT$_{Ox,fix}$=%6.1f" % (max_C, Tfuel_fixed, Toxid_fixed)
                    txt = plt.text(0.15, 0.60, txtstr, fontsize = 14, transform = fig.transFigure, verticalalignment = 'top', bbox = props)
                    plt.semilogx(a_max, T_max, 'b', figure = fig)
                    txt.set_visible(True)
                    plt.pause(0.05)
                    txt.set_visible(False)
                    fig1 = plt.figure(1)
                    plt.title(title_text, fontsize=10)
                    plt.plot(mixture_fraction_array, f.T, figure = fig1)

                # check for maximum number of flames
                n += 1
                if n > maxNumberOfFlames:
                    continuation = False
                    retrySolve = False;
                    break   
            
                # update step control 
                if stepSizeControl == True and flameControl == True:
                    stepFailed, dT = stepControl(ieval, N_opt, dT, MAXdT, MINdT)

            else:
                # failed to solve, prompt user for parameter modification or exit
                dlg = RetryDlg(OnePointControl, TwoPointControl, UserMaxTempPercentage, FuelSign, OxidSign, \
                               MINdT, MAXdT, dT, oxidizer_flux, cmean_skip)
                retrySolve = dlg.WaitResult()
                OnePointControl = dlg.one_point_control
                TwoPointControl = dlg.two_point_control
                UserMaxTempPercentage = dlg.max_temp
                FuelSign = dlg.fuel_sign
                OxidSign = dlg.oxidizer_sign
                MINdT = dlg.mindt
                MAXdT = dlg.maxdt
                dT = dlg.dt
                oxidizer_flux = dlg.mdot_oxidizer
                cmean_skip = dlg.cmean_skip
                restart = True

    # Calculating C = 0.0
    txt.set_visible(False)
    f.fuel_inlet.mdot *= 100
    f.oxidizer_inlet.mdot *=100
    f.set_max_grid_points(1, 101)
    for j in range(len(f.grid)):
        Tinter = min(Tfuel,Toxidizer)
        f.set_value(1, 'T', j, Tinter)	
    f.set_flame_control(1, StrainEq, UnityLewisNumber, False, False, -2000, -20, -2000, -20, False)
    f.solve(loglevel=loglevel, refine_grid=False, auto=False)
    T_max.append(max(f.T))
    a_max.append(f.strain_rate('max'))
    m_dots.append(f.fuel_inlet.mdot)
    C_array = get_C()
    max_C = np.abs(np.max(C_array))
    mixture_fraction_array = get_mixture_fraction()
    FastLaminarFlamelets(max_C, mixture_fraction_array)
    PrintFlush("%6d %10.1f %10.1f %10.2f %10s %10s %12d" \
            % (n, max(f.T), f.strain_rate('max'), dT, "lower", "False", f.flame.n_points))
    fig = plt.figure(0)
    title_text = '{}-{} T{}={:6.1f}K, T{}={:6.1f}K, P{}={:6.0f}Pa'.format(fuel, oxidizer, "$_{Fuel}$", Tfuel, "$_{Oxidizer}$", Toxidizer, "$_{Op}$",pressure)
    txtstr =  "$\Lambda$ = %6.4f\nT$_{F,fix}$=%6.1f\nT$_{Ox,fix}$=%6.1f" % (0.0, Tfuel, Toxidizer)
    txt = plt.text(0.15, 0.60, txtstr, fontsize = 14, transform = fig.transFigure, verticalalignment = 'top', bbox = props)
    plt.semilogx(a_max[:-1], T_max[:-1], 'b', figure = fig)
    plt.semilogx([a_max[-1], max(a_max)], [T_max[-1], T_max[-1]], 'b', figure = fig)
    txt.set_visible(True)
#    plt.pause(1000)
    fig1 = plt.figure(1)
    plt.plot(mixture_fraction_array, f.T, figure = fig1)

    # all done, close the file
    fast_file.close()
    summary.close()
    
    dlg = MessageDlg("1D Counterflow Flames Complete!")
    dlg.WaitResult()
    plt.close('all')


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
        PrintFlush('Run "PDFScript.py --help" to see options.')
        sys.exit(1)

    plt.ion()
    DiffusionFlame(options)
    sys.exit()

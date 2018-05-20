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

import os
import tkinter
import sys

from Utils import PrintFlush

###############################################################################
#
#  Defines class used to read the table definition
#
###############################################################################

class TableDef():

    def __init__(self):
        
        # initialize all member variables organized by tab
        
        # models tab
        self._reaction_mechanism = 'gri30.xml'
        self._model = "FPI"
        self._flow = "Turbulent"

        # boundary tab
        self._pressure = 101325.0 # pressure [Pa], information passed from C# dialog GUI
        self._Tfuel =  300.0 # unburned gas temperature [K], information pass from C# dialog GUI
        self._Toxidizer = 300.0 # unburned gas temperature [K], information pass from C# dialog GUI
        self._mole_fractions = True # this information is passed from C# dialog GUI
        self._fuel = ['CH4']  #the user chooses the fuel, this should be an array; this are the species which are greater than zero
                        # in the GUI
        self._fuel_fraction = [1.0] # this fraction is passed from the GUI
        self._oxidizer = ['O2', 'N2'] #the user chooses the oxidizer, this should be array; these are species > 0
        self._oxidizer_fraction = [0.21, 0.79] #these are the mole fractions chosen by user
        self._init_fuel_flux = 0.1 #kg/m^2-s, this is only a FPV input

        # numerics tab
        self._start = 0.0 # starting point (X0)
        self._end = 0.2 # ending point (X1)
        self._num_points = 21 # number of points between X0 and X1
        self._ss_rel = 1.0e-7 # steady state relative
        self._ss_abs = 1.0e-11 # steady state absolute
        self._tr_rel = 1.0e-7 # transient relative
        self._tr_abs = 1.0e-11 # transient absolute
        self._ss_age = 5  #Set the maximum number of times the Jacobian will be used before it must be re-evaluated by the steady solver.
        self._ts_age = 5  #Set the maximum number of times the Jacobian will be used before it must be re-evaluated by the transient solver. 
        self._stepsize = 1e-5 # initial time step in seconds
        self._nstep = [2,5,10,20] # sequence of integer numbers
        self._loglevel = 0 # amount of diagnostic output (0 to 8)
        self._refine = True  # 'True' to enable refinement, 'False' to disable
        self._temperature_limit_extinction =  2000 # K  #lets add this parameter to the numerics tab
        self._ratio = 10.0
        self._slope = 0.2
        self._curve = 0.2
        self._prune = 0.05
        self._expected_lean_extinction = 0.6 # fuel-lean extinction limit; THIS IS NOT NEEDED FOR FPV
        self._expected_rich_extinction = 6.0 # fuel-rich extinction limit; THIS IS NOT NEEDED FOR FPV

        # Table Generator tab
        self._composition=['CO','CO2','H2','H2O']  # progress variable definition
        self._z_points = 21 #number of fmean. This is the same as 21 in inital_grid variable above
        self._cmean_skip = 5 # number of cmean_skip for cmean entry in the GUI
        self._transport_model = 'Mix' # this can be either mix or multi or unity
        self._hmean = 0 # enthalpy levels
        self._nZVar = 0 # mixture fraction variance (fvar)
        self._nCVar = 0 # progress variance (cvar)
        self._z_component = ['H', 'O'] # list of elements
        self._pdf_type = "Beta"

        # Flamelet Control
        self._strain_rate = False
        self._precalculations = 1
        self._user_strain_factor = 1.0
        self._zero_order_steps = 500
        self._user_max_temp_percentage = 0.9
        self._fuel_sign = -1
        self._oxidizer_sign = -1
        self._dT = 20.0
        self._failed_solutions_to_terminate = 10
        self._max_num_flames = 10000
        self._nopt = 20
        self._max_dT = 20.0
        self._min_dT = 1.0
        self._one_point_control = False
        self._two_point_control = True
        self._oxidizer_mass_flux = 0.0

        # output tab
        self._YQOI = ['CO','CO2','H2','H2O','CH4','O2','N2','NO'] #user chooses species mass fraction to output
        self._RRQOI = ['CO','CO2','H2','H2O']  # user chooses species net reaction rates (kg/m^3-s) to output
        self.interpolation_type = 'slinear'


    @property
    def reaction_mechanism(self):
        '''reaction_mechanism'''
        return self._reaction_mechanism

    @property
    def model(self):
        '''model'''
        return self._model

    @property
    def flow(self):
        '''flow'''
        return self._flow

    @property
    def pressure(self):
        '''pressure'''
        return self._pressure

    @property
    def Tfuel(self):
        '''fuel temperature'''
        return self._Tfuel

    @property
    def Toxidizer(self):
        '''oxidizer temperature'''
        return self._Toxidizer

    @property
    def mole_fractions(self):
        '''mole fraction'''
        return self._mole_fractions

    @property
    def fuel(self):
        '''fuels array'''
        return self._fuel

    @property
    def fuel_fraction(self):
        '''fuel fraction array'''
        return self._fuel_fraction

    @property
    def oxidizer(self):
        '''oxidizers array'''
        return self._oxidizer

    @property
    def oxidizer_fraction(self):
        '''oxidizer fractions array'''
        return self._oxidizer_fraction

    @property
    def init_fuel_flux(self):
        '''initial fuel flux'''
        return self._init_fuel_flux
 
    @property
    def start(self):
        '''start (x0)'''
        return self._start
 
    @property
    def end(self):
        '''end (x1)'''
        return self._end
 
    @property
    def num_points(self):
        '''num_points'''
        return self._num_points
 
    @property
    def ss_rel(self):
        '''steady state relative (ss_rel)'''
        return self._ss_rel
 
    @property
    def ss_abs(self):
        '''steady state absolute (ss_abs)'''
        return self._ss_abs
 
    @property
    def tr_rel(self):
        '''transient relative (tr_rel)'''
        return self._tr_rel
 
    @property
    def tr_abs(self):
        '''transient absolute (tr_abs)'''
        return self._tr_abs
 
    @property
    def ss_age(self):
        '''steady state jacobian age (ss_age)'''
        return self._ss_age
 
    @property
    def ts_age(self):
        '''transient jacobian age (ts_age)'''
        return self._ts_age
 
    @property
    def stepsize(self):
        '''time step size'''
        return self._stepsize

    @property
    def nstep(self):
        '''time step nums (nstep)'''
        return self._nstep
 
    @property
    def loglevel(self):
        '''log level'''
        return self._loglevel
 
    @property
    def refine(self):
        '''refine grid'''
        return self._refine
 
    @property
    def temperature_limit_extinction(self):
        '''temperature limit extinction'''
        return self._temperature_limit_extinction
 
    @property
    def ratio(self):
        '''grid point spacing ratio (ratio)'''
        return self._ratio
 
    @property
    def slope(self):
        '''grid point spacing slope (slope)'''
        return self._slope
 
    @property
    def curve(self):
        '''grid point spacing curve (curve)'''
        return self._curve
 
    @property
    def prune(self):
        '''grid point spacing prune (prune)'''
        return self._prune
  
    @property
    def expected_lean_extinction(self):
        '''expected lean extinction'''
        return self._expected_lean_extinction
  
    @property
    def expected_rich_extinction(self):
        '''expected rich extinction'''
        return self._expected_rich_extinction
  
    @property
    def composition(self):
        '''progress variable definition array (composition)'''
        return self._composition
  
    @property
    def z_points(self):
        '''fmean (z_points)'''
        return self._z_points
  
    @property
    def cmean_skip(self):
        '''number progress points (cmean_skip)'''
        return self._cmean_skip
  
    @property
    def transport_model(self):
        '''transport model'''
        return self._transport_model
  
    @property
    def hmean(self):
        '''number enthalpy levels (hmean)'''
        return self._hmean
  
    @property
    def nZVar(self):
        '''mixture fraction variance (nZVar))'''
        return self._nZVar
   
    @property
    def nCVar(self):
        '''progress variance (nCVar)'''
        return self._nCVar
   
    @property
    def z_component(self):
        '''progress element list (z_component)'''
        return self._z_component

    @property
    def pdf_type(self):
        '''pdf distribution type'''
        return self._pdf_type
      
    @property
    def YQOI(self):
        '''array of species mass fraction to output (YQOI)'''
        return self._YQOI
    
    @property
    def RRQOI(self):
        '''array of species net reaction rates to output (RRQOI)'''
        return self._RRQOI
    
    @property
    def strain_rate(self):
        '''use strain equations flag'''
        return self._strain_rate
    
    @property
    def precalculations(self):
        '''initialization pressure'''
        return self._precalculations
    
    @property
    def user_strain_factor(self):
        '''user strain factor'''
        return self._user_strain_factor
    
    @property
    def zero_order_steps(self):
        '''max number of scaled zero order continuation'''
        return self._zero_order_steps
    
    @property
    def user_max_temp_percentage(self):
        '''max temperature percentage'''
        return self._user_max_temp_percentage
    
    @property
    def fuel_sign(self):
        '''increase or decrease fuel side'''
        return self._fuel_sign
    
    @property
    def oxidizer_sign(self):
        '''increase or decrease oxidizer side'''
        return self._oxidizer_sign
    
    @property
    def dT(self):
        '''temperature reduction delta'''
        return self._dT
    
    @property
    def failed_solutions_to_terminate(self):
        '''number of solution iterations before failing'''
        return self._failed_solutions_to_terminate
    
    @property
    def max_num_flames(self):
        '''maximum number of flames'''
        return self._max_num_flames
    
    @property
    def nopt(self):
        '''optimum number of newton iterations'''
        return self._nopt
    
    @property
    def max_dT(self):
        '''maximum temperature reduction'''
        return self._max_dT
    
    @property
    def min_dT(self):
        '''minmum temperature reduction'''
        return self._min_dT
    
    @property
    def one_point_control(self):
        '''one point control'''
        return self._one_point_control
    
    @property
    def two_point_control(self):
        '''two point control'''
        return self._two_point_control
    
    @property
    def oxidizer_mass_flux(self):
        '''oxidizer mass flux'''
        return self._oxidizer_mass_flux

    # method that reads in TableDefinition.txt file and updates member variables
    def Read(self):
        # read in fuels definition file
        PrintFlush("*** reading input file - TableDefinition.txt ***")
        try:
            file3=open('TableDefinition.txt', 'r')

        except Exception:
            PrintFlush("*** Error reading input file - TableDefinition.txt")
            PrintFlush("*** Make sure Working Directory set to path where python scripts reside")
            return False

        #read cantera file name
        temp = file3.readline().split("=")
        self._reaction_mechanism = temp[1].strip()

        #read flow
        temp = file3.readline().split("=")
        self._flow = temp[1].strip()

        # read mole fraction flag
        temp = file3.readline().split("=")
        fraction = temp[1].strip()
        if fraction == "Mole":
                self._mole_fractions = True
        else:
                self._mole_fractions = False

        # read pressure
        temp = file3.readline().split("=")
        self._pressure = float(temp[1].strip())

        # read fuel Temperature
        temp = file3.readline().split("=")
        self._Tfuel = float(temp[1].strip())

        # read oxidizer Temperature
        temp = file3.readline().split("=")
        self._Toxidizer = float(temp[1].strip())

        #read fuels
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        self._fuel = temp.split()

        #read fuel fraction
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        fraction_array = temp.split()
        self._fuel_fraction =  [float(i) for i in fraction_array]

        #read oxidizers
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        self._oxidizer = temp.split()

        #read oxidizer fraction
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        oxidizer_array = temp.split()
        self._oxidizer_fraction =  [float(i) for i in oxidizer_array]

        #read start (x0)
        temp = file3.readline().split("=")
        self._start = float(temp[1].strip())

        #read end (x1)
        temp = file3.readline().split("=")
        self._end = float(temp[1].strip())

        #read number of points
        temp = file3.readline().split("=")
        self._num_points = int(temp[1].strip())

        #read steady state absolute
        temp = file3.readline().split("=")
        self._ss_abs = float(temp[1].strip())

        #read steady state relative
        temp = file3.readline().split("=")
        self._ss_rel = float(temp[1].strip())

        #read transient absolute
        temp = file3.readline().split("=")
        self._tr_abs = float(temp[1].strip())

        #read transient relative
        temp = file3.readline().split("=")
        self._tr_rel = float(temp[1].strip())

        #read log level
        temp = file3.readline().split("=")
        self._loglevel = int(temp[1].strip())

        #read refine grid flag
        temp = file3.readline().split("=")
        itemp = int(temp[1].strip())
        if itemp == 0:
            self._refine = False
        else:
            self._refine = True
            
        #read steady state jacobian age
        temp = file3.readline().split("=")
        self._ss_age = int(temp[1].strip())
            
        #read transient jacobian age
        temp = file3.readline().split("=")
        self._ts_age = int(temp[1].strip())

        #read time step size
        temp = file3.readline().split("=")
        self._stepsize = float(temp[1].strip())

        #read time step nums
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        step_nums_array = temp.split()
        self._nstep =  [int(i) for i in step_nums_array]

        #read temperature limit extinction
        temp = file3.readline().split("=")
        self._temperature_limit_extinction = float(temp[1].strip())

        #read expected lean extinction
        temp = file3.readline().split("=")
        self._expected_lean_extinction = float(temp[1].strip())

        #read expected rich extinction
        temp = file3.readline().split("=")
        self._expected_rich_extinction = float(temp[1].strip())
 
        #read grid point spacing ratio
        temp = file3.readline().split("=")
        self._ratio = float(temp[1].strip())

        #read grid point slope
        temp = file3.readline().split("=")
        self._slope = float(temp[1].strip())

        #read grid point curve
        temp = file3.readline().split("=")
        self._curve = float(temp[1].strip())

        #read grid point prune
        temp = file3.readline().split("=")
        self._prune = float(temp[1].strip())

        #read species mass fraction to output
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        self._YQOI = temp.split()

        #read species net reaction rates (kg/m^3-s) to output
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        self._RRQOI = temp.split()

        #read progress variable definition
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        self._composition = temp.split()

        #read fmean
        temp = file3.readline().split("=")
        self._z_points = int(temp[1].strip())

        #read transport model
        temp = file3.readline().split("=")
        self._transport_model = temp[1].strip()

        #read initial fuel flux (FPV only)
        temp = file3.readline().split("=")
        self._init_fuel_flux = float(temp[1].strip())

        #read number of progress points (cmean)
        temp = file3.readline().split("=")
        self._cmean_skip = int(temp[1].strip())

        #read number of enthalpy levels (hmean)
        temp = file3.readline().split("=")
        self._hmean = int(temp[1].strip())

        #read mixture fraction variance (fvar)
        temp = file3.readline().split("=")
        self._nZVar = int(temp[1].strip())

        #read progress variance (cvar)
        temp = file3.readline().split("=")
        self._nCVar = int(temp[1].strip())
        
        #read element list
        temp = file3.readline().split("=")
        temp = temp[1].strip()
        self._z_component = temp.split()

        #read pdf type
        temp = file3.readline().split("=")
        self._pdf_type = temp[1].strip()

        #read strain rate flag
        temp = file3.readline().split("=")
        itemp = int(temp[1].strip())
        if itemp == 0:
            self._strain_rate = False
        else:
            self._strain_rate = True

        #read precalculations
        temp = file3.readline().split("=")
        self._precalculations = int(temp[1].strip())

        #read user_strain_factor
        temp = file3.readline().split("=")
        self._user_strain_factor = float(temp[1].strip())

        #read zero_order_steps
        temp = file3.readline().split("=")
        self._zero_order_steps = int(temp[1].strip())

        #read user_max_temp_percentage
        temp = file3.readline().split("=")
        self._user_max_temp_percentage = float(temp[1].strip())/100.0

        #read fuel_sign
        temp = file3.readline().split("=")
        self._fuel_sign = int(temp[1].strip())

        #read oxidizer_sign
        temp = file3.readline().split("=")
        self._oxidizer_sign = int(temp[1].strip())

        #read dT
        temp = file3.readline().split("=")
        self._dT = float(temp[1].strip())

        #read failed_solutions_to_terminate
        temp = file3.readline().split("=")
        self._failed_solutions_to_terminate = int(temp[1].strip())

        #read max_num_flames
        temp = file3.readline().split("=")
        self._max_num_flames = int(temp[1].strip())

        #read nopt
        temp = file3.readline().split("=")
        self._nopt = int(temp[1].strip())

        #read max_dT
        temp = file3.readline().split("=")
        self._max_dT = float(temp[1].strip())

        #read min_dT
        temp = file3.readline().split("=")
        self._min_dT = float(temp[1].strip())

        #read one point control flag flag
        temp = file3.readline().split("=")
        itemp = int(temp[1].strip())
        if itemp == 0:
            self._one_point_control = False
        else:
            self._one_point_control = True

        #read two point control flag flag
        temp = file3.readline().split("=")
        itemp = int(temp[1].strip())
        if itemp == 0:
            self._two_point_control = False
        else:
            self._two_point_control = True

        #read oxidizer_mass_flux
        temp = file3.readline().split("=")
        self._oxidizer_mass_flux = float(temp[1].strip())

        #read model
        temp = file3.readline().split("=")
        self._model = temp[1].strip()

        PrintFlush("*** input file read complete ***")
        file3.close()
        return True

     # method that prints contents of member variables
    def Report(self):
        PrintFlush("*** Reporting contents of TableDef ***")
        PrintFlush(TableDef.reaction_mechanism.__doc__ + " = ",self.reaction_mechanism)
        PrintFlush(TableDef.flow.__doc__ + " = ",self.flow)
        PrintFlush(TableDef.mole_fractions.__doc__ + " = ", self.mole_fractions)
        PrintFlush(TableDef.pressure.__doc__ + " = ", str(self.pressure))
        PrintFlush(TableDef.Tfuel.__doc__ + " = ", str(self.Tfuel))
        PrintFlush(TableDef.Toxidizer.__doc__ + " = ", str(self.Toxidizer))
        PrintFlush(TableDef.fuel.__doc__ + " = ", self.fuel)
        PrintFlush(TableDef.fuel_fraction.__doc__ + " = ", [str(i) for i in self.fuel_fraction])
        PrintFlush(TableDef.oxidizer.__doc__ + " = ", self.oxidizer)
        PrintFlush(TableDef.fuel_fraction.__doc__ + " = ", [str(i) for i in self.fuel_fraction])
        PrintFlush(TableDef.start.__doc__ + " = ", str(self.start))
        PrintFlush(TableDef.end.__doc__ + " = ", str(self.end))
        PrintFlush(TableDef.num_points.__doc__ + " = ", str(self.num_points))
        PrintFlush(TableDef.ss_abs.__doc__ + " = ", str(self.ss_abs))
        PrintFlush(TableDef.ss_rel.__doc__ + " = ", str(self.ss_rel))
        PrintFlush(TableDef.tr_abs.__doc__ + " = ", str(self.tr_abs))
        PrintFlush(TableDef.tr_rel.__doc__ + " = ", str(self.tr_rel))
        PrintFlush(TableDef.loglevel.__doc__ + " = ", str(self.loglevel))
        PrintFlush(TableDef.refine.__doc__ + " = ", self.refine)
        PrintFlush(TableDef.ss_age.__doc__ + " = ", str(self.ss_age))
        PrintFlush(TableDef.ts_age.__doc__ + " = ", str(self.ts_age))
        PrintFlush(TableDef.stepsize.__doc__ + " = ", str(self.stepsize))
        PrintFlush(TableDef.nstep.__doc__ + " = ", [str(i) for i in self.nstep])
        PrintFlush(TableDef.temperature_limit_extinction.__doc__ + " = ", str(self.temperature_limit_extinction))
        PrintFlush(TableDef.expected_lean_extinction.__doc__ + " = ", str(self.expected_lean_extinction))
        PrintFlush(TableDef.expected_rich_extinction.__doc__ + " = ", str(self.expected_rich_extinction))
        PrintFlush(TableDef.ratio.__doc__ + " = ", str(self.ratio))
        PrintFlush(TableDef.slope.__doc__ + " = ", str(self.slope))
        PrintFlush(TableDef.curve.__doc__ + " = ", str(self.curve))
        PrintFlush(TableDef.prune.__doc__ + " = ", str(self.prune))
        PrintFlush(TableDef.YQOI.__doc__ + " = ", self.YQOI)
        PrintFlush(TableDef.RRQOI.__doc__ + " = ", self.RRQOI)
        PrintFlush(TableDef.composition.__doc__ + " = ", self.composition)
        PrintFlush(TableDef.z_points.__doc__ + " = ", str(self.z_points))
        PrintFlush(TableDef.transport_model.__doc__ + " = ", self.transport_model)
        PrintFlush(TableDef.init_fuel_flux.__doc__ + " = ", str(self.init_fuel_flux))
        PrintFlush(TableDef.cmean_skip.__doc__ + " = ", str(self.cmean_skip))
        PrintFlush(TableDef.hmean.__doc__ + " = ", str(self.hmean))
        PrintFlush(TableDef.nZVar.__doc__ + " = ", str(self.nZVar))
        PrintFlush(TableDef.nCVar.__doc__ + " = ", str(self.nCVar))
        PrintFlush(TableDef.z_component.__doc__ + " = ", self.z_component)
        PrintFlush(TableDef.pdf_type.__doc__ + " = ", self._pdf_type)
        PrintFlush(TableDef.strain_rate.__doc__ + " = ", str(self.strain_rate))
        PrintFlush(TableDef.precalculations.__doc__ + " = ", str(self.precalculations))
        PrintFlush(TableDef.user_strain_factor.__doc__ + " = ", str(self.user_strain_factor))
        PrintFlush(TableDef.zero_order_steps.__doc__ + " = ", str(self.zero_order_steps))
        PrintFlush(TableDef.user_max_temp_percentage.__doc__ + " = ", str(self.user_max_temp_percentage))
        PrintFlush(TableDef.fuel_sign.__doc__ + " = ", str(self.fuel_sign))
        PrintFlush(TableDef.oxidizer_sign.__doc__ + " = ", str(self.oxidizer_sign))
        PrintFlush(TableDef.dT.__doc__ + " = ", str(self.dT))
        PrintFlush(TableDef.failed_solutions_to_terminate.__doc__ + " = ", str(self.failed_solutions_to_terminate))
        PrintFlush(TableDef.max_num_flames.__doc__ + " = ", str(self.max_num_flames))
        PrintFlush(TableDef.nopt.__doc__ + " = ", str(self.nopt))
        PrintFlush(TableDef.max_dT.__doc__ + " = ", str(self.max_dT))
        PrintFlush(TableDef.min_dT.__doc__ + " = ", str(self.min_dT))
        PrintFlush(TableDef.one_point_control.__doc__ + " = ", str(self.one_point_control))
        PrintFlush(TableDef.two_point_control.__doc__ + " = ", str(self.two_point_control))
        PrintFlush(TableDef.oxidizer_mass_flux.__doc__ + " = ", str(self.oxidizer_mass_flux))
        PrintFlush(TableDef.model.__doc__ + " = ", str(self.model))
        PrintFlush("*** Report Complete ***")

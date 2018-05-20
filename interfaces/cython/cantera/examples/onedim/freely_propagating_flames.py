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
import warnings
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import getopt

from table_def import TableDef
from Utils import MessageDlg, PrintFlush

def FreelyPropagatingFlames(option):

    # Suppress matplotlib warnings
    warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

    PrintFlush("Start flamelet generation of FPI type")

    # Set up a data directory for later reference
    data_directory = r"C:\WinPython-64bit-3.5.3.0Qt5\notebooks/"
    if not os.path.exists(data_directory):
        os.makedirs(data_directory)

    # Create a summary file for later reference
    file=open('premixed_flamelets.fla','w')
    summary=open('Summary.txt','w')
    summary.write('phi       SL      Tmax\n')

    # IdealGasMix object used to compute mixture properties, set to the state of the 
    # upstream fuel-air mixture

    # read in table definition file
    table = TableDef()
    if table.Read() != True:
        return
    table.Report()

    """
    ****************** begining of inputs ***************************************
    """

    """
    The following entry is obtained from the C# dialog models tab
    """

    gas = ct.Solution(table.reaction_mechanism)

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
    oxidizer_fraction = table.oxidizer_fraction #these are the mole fractions chosesn by user
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
    refine_grid = table.refine  # 'True' to enable refinement, 'False' to disable
    temperature_limit_extinction =  table.temperature_limit_extinction # K  #lets add this parameter to the numerics tab
    ratio = table.ratio
    slope = table.slope
    curve = table.curve
    prune = table.prune
    expected_lean_extinction = table.expected_lean_extinction # fuel-lean extinction limit
    expected_rich_extinction = table.expected_rich_extinction # fuel-rich extinction limit
    """
    end of boundary information
    """

    """
    The values below are coming from the C# dialog flamelet tab.
    For now only one entry is being passed because we are developing the laminar flamelet
    """
    composition = table.composition  # progress variable definition
    z_points = table.z_points #number of fmean. This is the same as 21 in inital_grid variable above
    transport_model = table.transport_model # this can be either mix or multi
    """
    end of flamelet information
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
    Tin = 0.5 * (Tfuel +Toxidizer)
    C_atoms = 0.0
    H_atoms = 0.0
    # Computing the stoichiometric fuel/air ratio (molar_sfar)
    for molecules in fuel:
        for i, char in enumerate(molecules):
            if char == 'C':
                flag = True
                C_atoms = 1
                if i + 1  < len(molecules):
                    for atoms in gas.element_names:
                        if molecules[i+1] == atoms:
                            flag = False
                    if flag == True:
                        C_atoms = int(molecules[i+1])
            if char == 'H':
                flag = True
                H_atoms = 1
                if i + 1 < len(molecules):
                    for atoms in gas.element_names:
                        if molecules[i+1] == atoms:
                            flag = False
                    if flag == True:
                        H_atoms = int(molecules[i+1])

    molar_sfar =  C_atoms + H_atoms / 4

    # Calculate the mass-based stoichiometric fuel air ratio (sfar)
    num = 0
    for i, molecules in enumerate(fuel):
        if mole_fractions:
            num += gas.molecular_weights[gas.species_index(fuel[i])] * fuel_fraction[i]
        else:
            num += fuel_fraction[i] / gas.molecular_weights[gas.species_index(fuel[i])]
    if mole_fractions == False:
        num /= num
    # Let's assume oxidation must have some amount of O2
    for i, molecules in enumerate(oxidizer):
        if molecules == 'O2':
            O2_index = i
                                     
    den = 0
    for i, molecules in enumerate(oxidizer):
        if mole_fractions:
            den += gas.molecular_weights[gas.species_index(oxidizer[i])] * oxidizer_fraction[i] / oxidizer_fraction[O2_index]
        else:
            den += oxidizer_fraction[i] / gas.molecular_weights[gas.species_index(fuel[i])]
    if mole_fractions == False:
        den /= den

    sfar= num / (molar_sfar * den)
    PrintFlush("Stoichiometric fuel ratio: ",sfar)

    # discretization of the mixture fraction (fmean) space
    npts = int(z_points/3)
    np_rich = 2*npts
    z_list = []
    for z_marker in range(0,z_points):
        if z_marker <= npts + 1:
            z_list.append(z_marker*sfar/npts)

        elif z_marker <= np_rich :
            z_list.append(sfar + sfar*(z_marker-npts)/(np_rich-npts))

        else:
            z_list.append(2*sfar + (1-2*sfar)*(z_marker-np_rich)/(z_points-1-np_rich))
    PrintFlush("Number of cases:",len(z_list))

    reactants_stoic = ''
    for i, molecules in enumerate(fuel):
        reactants_stoic += fuel[i] + ':' + str(fuel_fraction[i]) + ','
    for i, molecules in enumerate(oxidizer):
        reactants_stoic += oxidizer[i]+ ':' + str(molar_sfar * oxidizer_fraction[i] / oxidizer_fraction[O2_index])
        if len(oxidizer) > 1 and i < len(oxidizer) - 1:
            reactants_stoic += ','
    gas.TPX = Tin, pressure, reactants_stoic

    #Flame object
    f = ct.FreeFlame(gas, initial_grid)
    f.set_max_jac_age(ss_age,ts_age)
    f.set_time_step(stepsize,nstep)
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)

    if transport_model == 'Unity':
        PrintFlush("Unity = true")
        UnityLewisNumber = True
    else:
        PrintFlush("Unity = false")
        UnityLewisNumber = False

    f.set_flame_control(1, False, UnityLewisNumber, False, False, -2000, -20, -2000, -20)
    if refine_grid == True:
        f.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)

    """
    ********************end of Cantera settings for flamelet calculations***************************************
    """
    flamespeed_matrix = []
    T_max_matrix = []
    adiabatic_matrix = []

#    plt.ion()
    init_plot = True
    for index, mixture_fraction in enumerate(z_list):
        # premixed gas composition to Cantera
        if mixture_fraction != 1.0 and mixture_fraction != 0.0:
            phi = mixture_fraction / (1.0-mixture_fraction) / sfar
            reactants = ''
            for i, molecules in enumerate(fuel):
                reactants += fuel[i] + ':' + str(fuel_fraction[i]) + ','
            for i, molecules in enumerate(oxidizer):
                reactants += oxidizer[i]+ ':' + str((molar_sfar / phi )* oxidizer_fraction[i] / oxidizer_fraction[O2_index])
                if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                    reactants += ','
        elif mixture_fraction == 1.0:
            phi = 100000
            reactants = ''
            for i, molecules in enumerate(fuel):
                reactants += fuel[i] + ':' + str(fuel_fraction[i]) + ','
            for i, molecules in enumerate(oxidizer):
                reactants += oxidizer[i]+ ':0.0' 
                if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                    reactants += ','
        elif mixture_fraction == 0.0:
            phi = 0.0
            reactants = ''
            for i, molecules in enumerate(fuel):
                reactants += fuel[i] + ':0.0,'
            for i, molecules in enumerate(oxidizer):
                reactants += oxidizer[i]+ ':' + str(oxidizer_fraction[i] / oxidizer_fraction[O2_index]) 
                if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                    reactants += ','
        f.inlet.X = reactants
       
        # Initialize flame object
        gas.TPX = Tin, pressure, reactants
        gas.equilibrate('hp')    
        adiabatic_matrix.append(gas.T)
        gas.TPX = Tin, pressure, reactants
       
        try:
            #if phi > expected_lean_extinction and phi < expected_rich_extinction: 
            f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=True)           
            if transport_model == 'Multi':
                f.transport_model = 'Multi'
                f.solve(loglevel, refine_grid=refine_grid,auto=False)
                gas.TPX = Tin, pressure, reactants
                gas.equilibrate('hp')
            if f.u[0] > 0.0 and mixture_fraction > 0.0 and mixture_fraction < 1.0:
                PrintFlush('Flamespeed = %6.4f m/s for phi = %4.2f' % (f.u[0], phi))
            else:
                PrintFlush("Freely propagating premixed flame solution replaced with equilibrium calculation")
                #spurious_solution = True
                gas.TPX = Tin, pressure, reactants
                f = ct.FreeFlame(gas, np.linspace(table.start, table.end, 2))
                gas.equilibrate('hp')
                f.set_value(1, 'T', 1, gas.T)
                for j in range(gas.n_species):
                    f.set_value(1, 5 + j, 1, gas.Y[j])

        except Exception as e:
                PrintFlush("Freely propagating premixed flame solution replaced with equilibrium calculation")
                #spurious_solution = True
                gas.TPX = Tin, pressure, reactants
                f = ct.FreeFlame(gas, np.linspace(table.start, table.end, 2))
                gas.equilibrate('hp')
                f.set_value(1, 'T', 1, gas.T)
                for j in range(gas.n_species):
                    f.set_value(1, 5 + j, 1, gas.Y[j])

        flamespeed_matrix.append(f.u[0])
        T_max_matrix.append(f.T[-1])
        fig = plt.figure(0)
        fig.suptitle("Hybrid Equilibrium/Premixed Flames Manifold")
        title_text = '{}-{} T{}={:6.1f}K, T{}={:6.1f}K, P{}={:6.0f}Pa'.format(fuel, oxidizer, "$_{Fuel}$", Tfuel, "$_{Oxidizer}$", Toxidizer, "$_{Op}$",pressure)
        plt.title(title_text, fontsize=10)
        plt.plot(z_list[:index+1],T_max_matrix,c='b')
        plt.xlabel("Mixture Fraction, Z")
        plt.ylabel("Maximum Temperature, T$_{max}$ [K]")
        plt.pause(0.05)
        if init_plot:
        # only position initially, allows user to move afterwards
            mgr = plt.get_current_fig_manager()
            mgr.window.wm_geometry("+50+400")
        fig1 = plt.figure(1)
        fig1.suptitle("Equilibrium Calculation and Premixed Flames")
        plt.title(title_text, fontsize=10)
        plt.ylabel("Temperature, T [K]")
        plt.xlabel("Progress Variable, C")
        C_variable = []
        for i in range(len(f.grid)):
            reaction_progress=0.0
            for name in composition:
                reaction_progress += f.Y[gas.species_index(name),i]
            C_variable.append(reaction_progress)
        if init_plot:
        # only position initially, allows user to move afterwards
            mgr = plt.get_current_fig_manager()
            mgr.window.wm_geometry("+750+400")
            init_plot = False
        plt.plot(C_variable, f.T, figure = fig1)

        """
        ***********************Writes the flamelet table *****************************
        ********************************************************************************
        """   
        # Storing flamelet information
        file.write('HEADER\n')
        file.write('PREMIX_CHI 0.0000\n')
        #Proceeding with the output
        file.write('Z  %6.6f\n' % mixture_fraction)
        file.write('NUMOFSPECIES    %3d\n' % len(gas.Y))
        file.write('GRIDPOINTS    %3d\n' % (len(f.grid)))
        file.write('PRESSURE %6.2f\n' % gas.P)
        file.write('BODY\n')
        file.write('REACTION_PROGRESS\n')    
        # Calculate and write the progress variable
        C_variable = []
        for i in range(0,len(f.grid)):
            reaction_progress=0.0
            for name in composition:
                reaction_progress += f.Y[gas.species_index(name),i]
            C_variable.append(reaction_progress)
        max_C = max(C_variable)
    
        Lambda_variable = np.zeros(len(f.grid))
        for i in range(len(f.grid)):
            if i == 0:
                Lambda_variable[i] = 0.0
            elif i == len(f.grid)-1:
                Lambda_variable[i] = 1.0
            else:
                Lambda_variable[i] = C_variable[i]/max_C

        for i in range(0,len(f.grid)):
            file.write('%15.9e '% C_variable[i])
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                file.write('\n')   

        # Write temperature
        file.write('\nTEMPERATURE\n')
        for i in range(0,len(f.grid)):
            file.write('%15.9e '% f.T[i])
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                file.write('\n')            
   
        # Write heat capacity 
        file.write('\nCp\n')
        for i in range(0,len(f.grid)):
            file.write('%15.9e '% f.cp_mass[i])
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                file.write('\n')    
        
        # Write mean molecular weight 
        file.write('\nMMW\n')
        for i in range(0,len(f.grid)):
            file.write('%15.9e '% np.dot(f.X[:,i],gas.molecular_weights[:]))
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                file.write('\n')    

        # Write density 
        file.write('\nrho\n')
        for i in range(0,len(f.grid)):
            file.write('%15.9e '% f.density_mass[i])
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                file.write('\n')
        
        # Write viscosity 
        file.write('\nmu\n')
        for i in range(0,len(f.grid)):
            file.write('%15.9e '% f.viscosity[i])
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                file.write('\n')
        
        # Write thermal conductivity
        file.write('\nTC\n')
        for i in range(0,len(f.grid)):
            file.write('%15.9e '% f.thermal_conductivity[i])
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                file.write('\n')

        # Writing reaction progress source
        file.write('\nPREMIX_CDOT\n')
        for i in range(0,len(f.grid)):
            source=0.0
            for name in composition:
                source += (f.net_production_rates[gas.species_index(name),i] * gas.molecular_weights[gas.species_index(name)]) 
            file.write('%15.9e '% source)
            if ((i+1)%5==0 and (i+1) != len(f.grid)):
                    file.write('\n')                     
                
        # Write net production rates
        for name in RRQOI:
            file.write('\nnet_rxn-'+name+'\n')
            for i in range(0,len(f.grid)):
                source = f.net_production_rates[gas.species_index(name),i] * gas.molecular_weights[gas.species_index(name)]
                file.write('%15.9e '% source)
                if ((i+1)%5==0 and (i+1) != len(f.grid)):
                    file.write('\n')                

        # Write species mass fractions
        for name in YQOI:    
            file.write('\nmassfraction-'+name+'\n')
            for i in range(0,len(f.grid)):
                if f.Y[gas.species_index(name),i] < 1e-10 or ((name not in YQOI) and i == 0):  
                    file.write('%15.9e '% 0)
                else:
                    file.write('%15.9e '% f.Y[gas.species_index(name),i])
                if ((i+1)%5==0 and (i+1) != len(f.grid)):
                    file.write('\n')
                
        file.write('\n')
        file.write('\n')
    
        # Write a calculation summary
        summary.write('%6.3f  %6.3f    %6.2f\n' % (phi, f.u[0], np.max(f.T)))

        """
        ***********************the flamelet table completed ******************************
        ********************************************************************************
        """
    # Politely close all relevant files    
    file.close()
    summary.close()    

    dlg = MessageDlg("Freely Propagating Premixed Flames Complete!")
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
    FreelyPropagatingFlames(options)
    sys.exit()

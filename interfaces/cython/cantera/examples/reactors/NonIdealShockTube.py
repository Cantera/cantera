
# coding: utf-8

# # Non-Ideal Shock Tube Example
#  Ignition delay time computations in a high-pressure reflected shock tube reactor
# 
# In this example we illustrate how to setup and use a constant volume, adiabatic reactor to simulate
#   reflected shock tube experiments. This reactor will then be used to compute the ignition delay of
#   a gas at a specified initial temperature and pressure.  The example is written in a general way,
#   i.e., no particular EoS is presumed and ideal and real gas EoS can be used equally easily.
# 
# The reactor (system) is simply an 'insulated box,' and can technically be used for any number of
#   equations of state and constant-volume, adiabatic reactors.
# 
# Other than the typical Cantera dependencies, plotting functions require that you have matplotlib
#   installed, and data storing and analysis requires pandas. See https://matplotlib.org/ and
#   http://pandas.pydata.org/index.html, respectively, for additional info.

from __future__ import division
from __future__ import print_function

# Dependencies: pandas, numpy, and matplotlib.pyplot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager

import time

import cantera as ct
print('Runnning Cantera version: ' + ct.__version__)



# Define the ignition delay time (IDT).  This function computes the ignition delay from the occurence
#   of the peak concentration for the specified species.
def ignitionDelay(df, species):
    return df[species].argmax()


# Define the reactor temperature and pressure:
reactorTemperature = 1000 #Kelvin
reactorPressure = 40.0*101325.0 #Pascals

# Define the gas
# In this example we will choose a stoichiometric mixture of n-dodecane and air as the gas. For a
#   representative kinetic model, we use that developed by Wang, Ra, Jia, and Reitz
#   (https://www.erc.wisc.edu/chem_mech/nC12-PAH_mech.zip) by [H.Wang, Y.Ra, M.Jia, R.Reitz,
#   Development of a reduced n-dodecane-PAH mechanism. and its application for n-dodecane soot
#   predictions., Fuel 136 (2014) 25–36]

# R-K constants are calculated according to their critical temperature (Tc) and pressure (Pc):
#
#       a = 0.4275*(R^2)*(Tc^2.5)/(Pc)
#
#    and
#
#       b = 0.08664*R*Tc/Pc
#
#   where R is the gas constant.
#
#   For stable species, the critical properties are readily available.  For radicals and other
#   short-lived intermediates, the Joback method is used to estimate critical properties.  See Joback
#   and Reid, "Estimation of pur-component properties from group-contributions," Chem. Eng. Comm. 57
#   (1987) 233-243, for details of the method.

# There is a slight discontinuity in the thermo for three species at the mid-point temperatrue.  We
#   are aware and okay, so we will suppress the warning statement (note: use this feature at your own
#   risk, in other codes!)
ct.suppress_thermo_warnings()


"""Real gas IDT calculation"""
# Load the real gas mechanism:
real_gas = ct.Solution('nDodecane_Reitz.cti','nDodecane_RK')

# Set the state of the gas object:
real_gas.TP = reactorTemperature, reactorPressure

# Define the fuel, oxidizer and set the stoichiometry:
real_gas.set_equivalence_ratio(phi=1.0, fuel='c12h26',
                               oxidizer={'o2':1.0, 'n2':3.76})

# Create a reactor object and add it to a reactor network
# In this example, this will be the only reactor in the network
r = ct.Reactor(contents=real_gas)
reactorNetwork = ct.ReactorNet([r])

# now compile a list of all variables for which we will store data
stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]

# Use the above list to create a DataFrame
timeHistory_RG = pd.DataFrame(columns=stateVariableNames)
#Tic
t0 = time.time()

# This is a starting estimate. If you do not get an ignition within this time, increase it
estimatedIgnitionDelayTime = 0.005
t = 0

counter = 1;
while(t < estimatedIgnitionDelayTime):
    t = reactorNetwork.step()
    if (counter%20 == 0):
        # We will save only every 20th value. Otherwise, this takes too long
        # Note that the species concentrations are mass fractions
        timeHistory_RG.loc[t] = reactorNetwork.get_state()
    counter+=1

# We will use the 'oh' species to compute the ignition delay
tau_RG = ignitionDelay(timeHistory_RG, 'oh')

#Toc
t1 = time.time()
print('Computed Real Gas Ignition Delay: {:.3e} seconds. Took {:3.2f}s to compute'.format(tau_RG, t1-t0))


"""Ideal gas IDT calculation"""
# Create the ideal gas object:
ideal_gas = ct.Solution('nDodecane_Reitz.cti','nDodecane_IG')

# Set the state of the gas object:
ideal_gas.TP = reactorTemperature, reactorPressure

# Define the fuel, oxidizer and set the stoichiometry:
ideal_gas.set_equivalence_ratio(phi=1.0, fuel='c12h26',
                                oxidizer={'o2':1.0, 'n2':3.76})

r = ct.Reactor(contents=ideal_gas)
reactorNetwork = ct.ReactorNet([r])

# now compile a list of all variables for which we will store data
stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]

# Use the above list to create a DataFrame
timeHistory_IG = pd.DataFrame(columns=stateVariableNames)

#Tic
t0 = time.time()

t = 0

counter = 1;
while(t < estimatedIgnitionDelayTime):
    t = reactorNetwork.step()
    if (counter%20 == 0):
        # We will save only every 20th value. Otherwise, this takes too long
        # Note that the species concentrations are mass fractions
        timeHistory_IG.loc[t] = reactorNetwork.get_state()
    counter+=1

# We will use the 'oh' species to compute the ignition delay
tau_IG = ignitionDelay(timeHistory_IG, 'oh')

#Toc
t1 = time.time()

print('Computed Ideal Gas Ignition Delay: {:.3e} seconds. Took {:3.2f}s to compute'.format(tau_IG, t1-t0))
print('Ideal gas error: {:2.2f} %'.format(100*(tau_IG-tau_RG)/tau_RG))

# If you want to save all the data - molefractions, temperature, pressure, etc
# uncomment the next line
# timeHistory.to_csv("time_history.csv")

# Plot the result

plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['figure.autolayout'] = True

# Figure illustrating the definition of ignition delay time (IDT).

plt.figure()
plt.plot(timeHistory_RG.index, timeHistory_RG['oh'],'-o',color='b',markersize=4)
plt.plot(timeHistory_IG.index, timeHistory_IG['oh'],'-o',color='r',markersize=4)
plt.xlabel('Time (s)',fontname='Times New Roman')
plt.ylabel('$\mathdefault{OH\, mass\, fraction,}\,  \mathdefault{y_{OH}}$',
           fontname='Times New Roman')

# Figure formatting:
plt.xlim([0,0.00055])
ax = plt.gca()
font = plt.matplotlib.font_manager.FontProperties(family='Times New Roman',size=14)
ax.annotate("",xy=(tau_RG,0.005), xytext=(0,0.005),
            arrowprops=dict(arrowstyle="<|-|>",color='r',linewidth=2.0),
            fontsize=14,)
plt.annotate('Ignition Delay Time (IDT)', xy=(0,0), xytext=(0.00008, 0.00525),
             family='Times New Roman',fontsize=16);

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(12)
    tick.label1.set_fontname('Times New Roman')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(12)
    tick.label1.set_fontname('Times New Roman')

plt.legend(['Real Gas','Ideal Gas'],prop=font,frameon=0)

# If you want to save the plot, uncomment this line (and edit as you see fit):
#plt.savefig('IDT_nDodecane_1000K_40atm.pdf',dpi=350,format='pdf')


"""Demonstration of NTC behavior"""
# Let us use the reactor model to demonstrate the impacts of non-ideal behavior on IDTs in the
#   Negative Temperature Coefficient (NTC) region, where observed IDTs, counter to intuition, increase
#   with increasing temperature.

# Make a list of all the temperatures at which we would like to run simulations:
T = [1250, 1225, 1200, 1150, 1100, 1075, 1050, 1025, 1012.5, 1000, 987.5, 975, 962.5, 950,
     937.5, 925, 912.5, 900, 875, 850, 825, 800]

# If we desire, we can define different IDT starting guesses for each temperature:
estimatedIgnitionDelayTimes = np.ones(len(T))
# But we won't, at least in this example :)
estimatedIgnitionDelayTimes[:] = 0.005

# Now create a dataFrame for the real gas results:
ignitionDelays_RG = pd.DataFrame(data={'T':T})
ignitionDelays_RG['ignDelay'] = np.nan


# Now, we simply run the code above for each temperature.
"""Real Gas"""
for i, temperature in enumerate(T):
    # Setup the gas and reactor
    reactorTemperature = temperature
    real_gas.TP = reactorTemperature, reactorPressure
    real_gas.set_equivalence_ratio(phi=1.0, fuel='c12h26',
                                   oxidizer={'o2':1.0, 'n2':3.76})
    r = ct.Reactor(contents=real_gas)
    reactorNetwork = ct.ReactorNet([r])

    # Create and empty data frame
    timeHistory = pd.DataFrame(columns=timeHistory_RG.columns)

    t0 = time.time()

    t = 0
    counter = 0
    while t < estimatedIgnitionDelayTimes[i]:
        t = reactorNetwork.step()
        if not counter % 20:
            timeHistory.loc[t] = r.get_state()
        counter += 1

    tau = ignitionDelay(timeHistory, 'oh')
    t1 = time.time()

    print('Computed Real Gas Ignition Delay: {:.3e} seconds for T={}K. Took {:3.2f}s to compute'.format(tau, temperature, t1-t0))

    ignitionDelays_RG.set_value(index=i, col='ignDelay', value=tau)



"""Repeat for Ideal Gas"""
# Create a dataFrame for the ideal gas results:
ignitionDelays_IG = pd.DataFrame(data={'T':T})
ignitionDelays_IG['ignDelay'] = np.nan

for i, temperature in enumerate(T):
    # Setup the gas and reactor
    reactorTemperature = temperature
    ideal_gas.TP = reactorTemperature, reactorPressure
    ideal_gas.set_equivalence_ratio(phi=1.0, fuel='c12h26',
                                    oxidizer={'o2':1.0, 'n2':3.76})
    r = ct.Reactor(contents=ideal_gas)
    reactorNetwork = ct.ReactorNet([r])

    # Create and empty data frame
    timeHistory = pd.DataFrame(columns=timeHistory_IG.columns)

    t0 = time.time()

    t = 0
    counter = 0
    while t < estimatedIgnitionDelayTimes[i]:
        t = reactorNetwork.step()
        if not counter % 20:
            timeHistory.loc[t] = r.get_state()
        counter += 1

    tau = ignitionDelay(timeHistory, 'oh')
    t1 = time.time()

    print('Computed Ideal Gas Ignition Delay: {:.3e} seconds for T={}K. Took {:3.2f}s to compute'.format(tau, temperature, t1-t0))

    ignitionDelays_IG.set_value(index=i, col='ignDelay', value=tau)


# Figure: ignition delay ($\tau$) vs. the inverse of temperature ($\frac{1000}{T}$).
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(1000/ignitionDelays_RG['T'], 1e6*ignitionDelays_RG['ignDelay'],'-',
        linewidth=2.0,color='b')
ax.plot(1000/ignitionDelays_IG['T'], 1e6*ignitionDelays_IG['ignDelay'],'-.',
        linewidth=2.0,color='r')
ax.set_ylabel(r'$\mathdefault{Ignition\, Delay\, (\mu s)}$',fontname='Times New Roman',
              fontsize=16)
ax.set_xlabel(r'$\mathdefault{1000/T\, (K^{-1})}$',fontname='Times New Roman', fontsize=16)

ax.set_xlim([0.8,1.2])

# Add a second axis on top to plot the temperature for better readability
ax2 = ax.twiny()
ticks = ax.get_xticks()
ax2.set_xticks(ticks)
ax2.set_xticklabels((1000/ticks).round(1))
ax2.set_xlim(ax.get_xlim())
ax2.set_xlabel('Temperature (K)',fontname='Times New Roman',fontsize=16);

ticks_font = font_manager.FontProperties(family='Times New Roman', style='normal',
                                         size=12, weight='normal',
                                         stretch='normal')
for label in ax.get_yticklabels():
    label.set_fontproperties(ticks_font)
for label in ax.get_xticklabels():
    label.set_fontproperties(ticks_font)
for label in ax2.get_xticklabels():
    label.set_fontproperties(ticks_font)

ax.legend(['Real Gas','Ideal Gas'],prop=font,frameon=0,loc=2)

# If you want to save the plot, uncomment this line (and edit as you see fit):
#plt.savefig('NTC_nDodecane_40atm.pdf',dpi=350,format='pdf')

# Show the plots.
plt.show()





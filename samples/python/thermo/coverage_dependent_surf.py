"""
This example demonstrates (1). the four different dependency models available
for coverage-dependent enthalpy and entropy calculations and (2). capability
of including self-interaction but also cross-interaction among different
surface species.

First demonstration is with a binary system composed of Pt and CO* where Pt
is vacant Pt sites and CO* is CO adsorbates with CO*-CO* lateral interaction.
Coverage-dependent enthalpy and entropy calculated with the four dependency
models are plotted against discrete data points used for dependency model
parametrization.

Second demonstration is with a ternary system composed of Pt, CO* and O*.
In this system, enthalpy of CO* is dependent both of CO* coverage and O*
coverage. That is, it encompasses CO*-CO* interaction as well as CO*-O*
interaction. The CO* enthalpy is plotted as a function of CO* and O*
coverages.

Requires: cantera >= 3.0.0, matplotlib >= 2.0
Keywords: thermodynamics, surface chemistry, catalysis
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# first demonstration: the four coverage dependency models

# provide discrete enthalpy and entropy values calculated with DFT
# array of CO* coverage
dft_covs = [0., 0.11, 0.22, 0.33, 0.44, 0.56, 0.67, 0.78, 0.89, 1.0]
# array of nondimensionalized DFT-derived CO* enthalpy values
dft_hrts = [-112.8 , -112.8 , -108.78, -105.76, -106.11, -102.51,  -96.83,
            -85.03,  -78.76,  -72.19]
# array of nondimensionalized DFT-derived CO* entropy values
dft_srs = [2.46, 2.46, 3.79, 3.28, 3.14, 2.41, 2.19, 1.63, 1.6 , 1.12]

# create a dictionary to store coverage-dependent enthalpy and entropy
# from the four dependency models. key is the model name and value is
# a dictionary of yaml's phase name, enthalpy, and entropy
models_dict = {'linear': {'phase_name': 'covdep_lin'},
               'piecewise-linear': {'phase_name': 'covdep_pwlin'},
               'polynomial': {'phase_name': 'covdep_poly'},
               'interpolative': {'phase_name': 'covdep_int'}}

# create an array of CO* coverages from 0 (no CO*) to 1 (full mono-layer)
CO_covs = np.linspace(0, 1.0, 101)
# remember CO* species index
i_CO = 1

# get coverage-dependent enthalpy and entropy for all four models
for model in models_dict:
    # dependency parameters are provided in covdepsurf.yaml
    # import a phase from covdepsurf.yaml and set temperature and pressure
    phase = ct.Interface('covdepsurf.yaml',
                         models_dict[model]['phase_name'])
    phase.TP = 300, ct.one_atm
    # under the current model, create empty lists to store enthalpy and
    # entropy
    models_dict[model]['hrt'] = []
    models_dict[model]['sr'] = []
    # sweep through the array of CO* coverages and store enthalpy and entropy
    for CO_cov in CO_covs:
        # set coverages in order of [Pt_cov, CO_cov]
        phase.coverages = [1 - CO_cov, CO_cov]
        # get and store enthalpy and entropy values to the corresponding list
        models_dict[model]['hrt'].append(phase.standard_enthalpies_RT[i_CO])
        models_dict[model]['sr'].append(phase.standard_entropies_R[i_CO])

# plot coverage-dependent enthalpy against DFT-derived enthalpy data
plt.plot(dft_covs, dft_hrts, marker='o', color='k', linewidth=0,
         label='DFT data')
for model in models_dict:
    plt.plot(CO_covs, models_dict[model]['hrt'], label=model)
plt.xlabel(r'$\theta_\mathrm{CO*}$')
plt.ylabel(r'$h^{\circ}_\mathrm{CO*}/\mathrm{R}T$')
plt.legend(loc='best')
plt.title('CO* standard state enthalpy\n'
           + 'with four dependency models at 300 K, 1 atm')
plt.show()
plt.clf()

# plot coverage-dependent entropy against DFT-derived enthalpy data
plt.plot(dft_covs, dft_srs, marker='o', color='k', linewidth=0,
         label='DFT data')
for model in models_dict:
    plt.plot(CO_covs, models_dict[model]['sr'], label=model)
plt.xlabel(r'$\theta_\mathrm{CO*}$')
plt.ylabel(r'$s^{\circ}_\mathrm{CO*}/\mathrm{R}$')
plt.legend(loc='best')
plt.title('CO* standard state entropy\n'
          + 'with four dependency models at 300 K, 1 atm')
plt.show()
plt.clf()

# second demonstration: coverage dependence from self- and cross-interaction

# import a phase with both types of lateral interactions
phase = ct.Interface('covdepsurf.yaml', 'covdep_cross')
phase.TP = 300, ct.one_atm

# create an array of coverages
covs = np.linspace(0, 1.0, 101)
# create N by N matrix to store CO* enthalpy values
CO_hrts = np.zeros(shape=(len(covs),len(covs)))

# iterate through the coverage array and fill in CO*enthalpy values
for i, CO_cov in enumerate(covs):
    for j in range(len(covs)-i):
        O_cov = covs[j]
        # set coverages in order of [Pt_cov, CO_cov, O_cov]
        phase.coverages = [1. - CO_cov - O_cov, CO_cov, O_cov]
        CO_hrts[i,j] = phase.standard_enthalpies_RT[i_CO]

# plot the enthalpy matrix
cntr = plt.contourf(covs, covs, CO_hrts,
                    levels=np.linspace(-120, -50, 71),
                    cmap='inferno')

plt.xlabel(r'$\theta_\mathrm{CO*}$')
plt.ylabel(r'$\theta_\mathrm{O*}$')
cbar = plt.colorbar(cntr)
cbar.ax.set_ylabel(r'$h^{\circ}_\mathrm{CO*}/\mathrm{R}T$')
plt.title('CO* standard state enthalpy\n'
          + 'with CO*$-$CO* and CO*$-$O* interactions')
plt.show()

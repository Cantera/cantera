"""
One-dimensional packed-bed, catalytic-membrane reactor
======================================================

The model shown in this example simulates heterogeneous catalytic processes inside
packed-bed, catalytic membrane reactors. The gas-phase and surface-phase species
conservation equations are derived and the system of differential-algebraic equations
(DAE) is solved using the ``scikits.odes.dae`` IDA solver.

Requires: cantera >= 3.1.0, matplotlib >= 2.0, scikits.odes >= 2.7.0

.. tags:: Python, surface chemistry, 1D flow, catalysis, packed bed reactor,
          porous media, user-defined model
"""


# %%
# Methodology
# -----------
#
# A one-dimensional, steady-state catalytic-membrane reactor model with surface
# chemistry is developed to analyze species profiles along the length of a packed-bed,
# catalytic membrane reactor. The same model can further be simplified to simulate a
# simple packed-bed reactor by excluding the membrane. The example here demonstrates the
# one-dimensional reactor model explained by G. Kogekar [1]_.

# %%
# Governing equations
# ^^^^^^^^^^^^^^^^^^^
#
# Assuming steady-state, one-dimensional flow within the packed bed, total-mass, species
# mass and energy conservation may be stated as [2]_:
#
# .. math::
#
#     \frac{d(\rho u)}{dz} &= \sum_{k=1}^{K_\t{g}} \dot {s}_k W_k A_\t{s}
#                          + \frac{P_\t{b}}{A_\t{b}} j_{k_\t{M}}, \\
#
#     \rho u \frac{dY_k}{dz} + A_\t{s} Y_k \sum_{k=1}^{K_\t{g}} \dot {s}_k W_k &=
#          A_\t{s} \dot {s}_k W_k + \delta_{k, k_\t{M}} \frac{P_\t{b}}{A_\t{b}} j_{k_\t{M}}, \\
#
#     \rho u c_\t{p} \frac{dT}{dz}
#         + \sum_{k=1}^{K_\t{g}} h_k (\phi_\t{g} \dot {\omega}_k + A_\t{s} \dot {s}_k) W_k
#         &= \hat h \frac{P_\t{b}}{A_\t{b}}(T_\t{w} - T)
#         + \delta_{k, k_\t{M}} \frac{P_\t{b}}{A_\t{b}} h_{k_\t{M}}  j_{k_\t{M}}.
#
# The fractional coverages of the :math:`K_\t{s}` surface adsorbates :math:`\theta_k`
# must satisfy
#
# .. math::
#     \dot {s}_k = 0, \qquad (k = 1,\ldots, K_\t{s}),
#
# which, at steady state, requires no net production/consumption of surface species by
# the heterogeneous reactions [3]_.
#
# The pressure within the bed is calculated as:
#
# .. math::
#     \frac{dp}{dz} = - \left( \frac{\phi_\t{g} \mu}{\beta_\t{g}} \right) u,
#
# where :math:`\mu` is the gas-phase mixture viscosity. The packed-bed permeability
# :math:\beta_\t{g}` is evaluated using the Kozeny-Carman relationship as
#
# .. math::
#     \beta_\t{g} = \frac{\phi_\t{g}^3 D_\t{p}^2}{72 \tau_\t{g}(1 - \phi_\t{g})^2},
#
# where :math:`\phi_\t{g}`, :math:`\tau_\t{g}`, and :math:`D_\t{p}` are the bed
# porosity, tortuosity, and particle diameter, respectively.
#
# The independent variable in these conservation equations is the position :math:`z`
# along the reactor length. The dependent variables include total mass flux :math:`\rho
# u`, pressure :math:`p`, temperature :math:`T`, gas-phase mass fractions :math:`Y_k`,
# and surfaces coverages :math:`\theta_k`. Gas-phase fluxes through the membrane are
# represented as :math:`j_{k, \t{M}}`. Geometric parameters :math:`A_\t{s}`,
# :math:`P_\t{b}`, and :math:`A_\t{b}` represent the catalyst specific surface area
# (dimensions of surface area per unit volume), reactor perimeter, and reactor
# cross-sectional flow area, respectively. Other parameters include bed porosity
# :math:`\phi_\t{g}` and gas-phase species molecular weights :math:`W_k`. The gas
# density :math:`\rho` is evaluated using the equation of state (such as ideal gas,
# Redlich-Kwong or Peng-Robinson).
#
# If a perm-selective membrane is present, then :math:`j_{k_\t{M}}` represents the
# gas-phase flux through the membrane and :math:`{k_\t{M}}` is the gas-phase species
# that permeates through the membrane. The Kronecker delta, :math:`\delta_{k, k_\t{M}} =
# 1` for the membrane-permeable species and :math:`\delta_{k, k_\t{M}} = 0` otherwise.
# The membrane flux is calculated using Sievert's law as
#
# .. math::
#     j_{k_\t{M}}^{\text{Mem}} =
#         \frac{B_{k_\t{M}}}{t} \left( p_{k_\t{M} {\text{, mem}}}^\alpha - p_{k_\t{M} \text{, sweep}}^\alpha \right) W_{k_\t{M}}
#
# where :math:`B_{k_\t{M}}` is the membrane permeability, :math:`t` is the membrane
# thickness. :math:`p_{k_\t{M} \text{, mem}}` and :math:`p_{k_\t{M} \text{, sweep}}`
# represent perm-selective species partial pressures within the packed-bed and the
# exterior sweep channel. The present model takes the pressure exponent :math:`\alpha`
# to be unity. The membrane flux for all other species (:math:`k \neq k_\t{M}`) is zero.
#
# Chemistry mechanism
# ^^^^^^^^^^^^^^^^^^^
#
# This example uses a detailed 12-step elementary micro-kinetic reaction mechanism that
# describes ammonia formation and decomposition kinetics over the Ru/Ba-YSZ catalyst.
# The reaction mechanism is developed and validated using measured performance in a
# laboratory-scale packed-bed reactor [4]_. This example also incorporates a Pd-based
# H₂ perm-selective membrane.
#
# Solver
# ^^^^^^
#
# The above governing equations represent a complete solution for a steady-state
# packed-bed, catalytic membrane reactor model. The dependent variables are the
# mass-flux :math:`\rho u`, species mass-fractions :math:`Y_k`, pressure :math:`p`,
# temperature :math:`T`, and surface coverages :math:`\theta_k`. The equation of state
# is used to obtain the mass density, :math:`\rho`.
#
# The governing equations form an initial value problem (IVP) in a
# differential-algebraic (DAE) form as follows:
#
# .. math::
#    f(z,{\bf{y}}, {\bf y'}, {\bf c}) = {\bf 0},
#
# where :math:`\bf y` and :math:`\bf y'` represent the state vector and its derivative
# vector, respectively. All other constants such as reference temperature, chemical
# constants, etc. are incorporated in vector :math:`\bf c` (Refer to [1]_ for more
# details). This type of DAE system in this example is solved using the
# ``scikits.odes.dae`` IDA solver.

# %%
# Import Cantera and scikits
# --------------------------

import numpy as np
from scikits.odes import dae
import cantera as ct
import matplotlib.pyplot as plt

# %%
# Define gas-phase and surface-phase species
# ------------------------------------------

# Import the reaction mechanism for Ammonia synthesis/decomposition on Ru-Ba/YSZ catalyst
mechfile = "example_data/ammonia-Ru-Ba-YSZ-CSM-2019.yaml"
# Import the models for surface-phase and gas
surf = ct.Interface(mechfile, "Ru_surface")
gas = surf.adjacent["gas"]

# Other parameters
n_gas = gas.n_species  # number of gas species
n_surf = surf.n_species  # number of surface species
n_gas_reactions = gas.n_reactions  # number of gas-phase reactions

# Set offsets of dependent variables in the solution vector
offset_rhou = 0
offset_p = 1
offset_T = 2
offset_Y = 3
offset_Z = offset_Y + n_gas
n_var = offset_Z + n_surf  # total number of variables (rhou, P, T, Yk and Zk)

print("Number of gas-phase species = ", n_gas)
print("Number of surface-phase species = ", n_surf)
print("Number of variables = ", n_var)

# %%
# Define reactor geometry and operating conditions
# ------------------------------------------------

# Reactor geometry
L = 5e-2  # length of the reactor (m)
R = 5e-3  # radius of the reactor channel (m)
phi = 0.5  # porosity of the bed (-)
tau = 2.0  # tortuosity of the bed (-)
dp = 3.37e-4  # particle diameter (m)
As = 3.5e6  # specific surface area (1/m)

# Energy (adiabatic or isothermal)
solve_energy = True  # True: Adiabatic, False: isothermal

# Membrane (True: membrane, False: no membrane)
membrane_present = True
membrane_perm = 1e-15  # membrane permeability (kmol*m3/s/Pa)
thickness = 3e-6  # membrane thickness (m)
membrane_sp_name = "H2"  # membrane-permeable species name
p_sweep = 1e5  # partial pressure of permeable species in the sweep channel (Pa)
permeance = membrane_perm / thickness  # permeance of the membrane (kmol*m2/s/Pa)

if membrane_present:
    print("Modeling packed-bed, catalytic-membrane reactor...")
    print(membrane_sp_name, "permeable membrane is present.")

# Get required properties based on the geometry and mechanism
W_g = gas.molecular_weights  # vector of molecular weight of gas species
vol_ratio = phi / (1 - phi)
eff_factor = phi / tau  # effective factor for permeability calculation
# permeability based on Kozeny-Carman equation
B_g = B_g = vol_ratio**2 * dp**2 * eff_factor / 72
area2vol = 2 / R  # area to volume ratio assuming a cylindrical reactor
D_h = 2 * R  # hydraulic diameter
membrane_sp_ind = gas.species_index(membrane_sp_name)

# Inlet operating conditions
T_in = 673  # inlet temperature [K]
p_in = 5e5  # inlet pressure [Pa]
v_in = 0.001  # inlet velocity [m/s]
T_wall = 723  # wall temperature [K]
h_coeff = 1e2  # heat transfer coefficient [W/m2/K]

# Set gas and surface states
gas.TPX = T_in, p_in, "NH3:0.99, AR:0.01"  # inlet composition
surf.TP = T_in, p_in
Yk_0 = gas.Y
rhou0 = gas.density * v_in

# Initial surface coverages
# advancing coverages over a long period of time to get the steady state.
surf.advance_coverages(1e10)
Zk_0 = surf.coverages


# %%
# Define residual function required for IDA solver
# ------------------------------------------------

def residual(z, y, yPrime, res):
    """Solution vector for the model
    y = [rho*u, p, T, Yk, Zk]
    yPrime = [d(rho*u)dz, dpdz, dTdz, dYkdz, dZkdz]
    """
    # Get current thermodynamic state from solution vector and save it to local variables.
    rhou = y[offset_rhou]  # mass flux (density * velocity)
    Y = y[offset_Y : offset_Y + n_gas]  # vector of mass fractions
    Z = y[offset_Z : offset_Z + n_surf] # vector of site fractions
    p = y[offset_p]  # pressure
    T = y[offset_T]  # temperature

    # Get derivatives of dependent variables
    drhoudz = yPrime[offset_rhou]
    dYdz = yPrime[offset_Y : offset_Y + n_gas]
    dZdz = yPrime[offset_Z : offset_Z + n_surf]
    dpdz = yPrime[offset_p]
    dTdz = yPrime[offset_T]

    # Set current thermodynamic state for the gas and surface phases
    # Note: use unnormalized mass fractions and site fractions to avoid
    # over-constraining the system
    gas.set_unnormalized_mass_fractions(Y)
    gas.TP = T, p
    surf.set_unnormalized_coverages(Z)
    surf.TP = T, p

    # Calculate required variables based on the current state
    coverages = surf.coverages  # surface site coverages
    # heterogeneous production rate of gas species
    sdot_g = surf.get_net_production_rates("gas")
    # heterogeneous production rate of surface species
    sdot_s = surf.get_net_production_rates("Ru_surface")
    wdot_g = np.zeros(n_gas)
    # specific heat of the mixture
    cp = gas.cp_mass
    # partial enthalpies of gas species
    hk_g = gas.partial_molar_enthalpies

    if n_gas_reactions > 0:
        # homogeneous production rate of gas species
        wdot_g = gas.net_production_rates
    mu = gas.viscosity  # viscosity of the gas-phase

    # Calculate density using equation of state
    rho = gas.density

    # Calculate flux term through the membrane
    # partial pressure of membrane-permeable species
    memsp_pres = p * gas.X[membrane_sp_ind]
    # negative sign indicates the flux going out
    membrane_flux = -permeance * (memsp_pres - p_sweep) * W_g[membrane_sp_ind]

    # Conservation of total-mass
    # temporary variable
    sum_continuity = As * np.sum(sdot_g * W_g) + phi * np.sum(wdot_g * W_g)
    res[offset_rhou] = (drhoudz - sum_continuity
                        - area2vol * membrane_flux * membrane_present)

    # Conservation of gas-phase species
    res[offset_Y:offset_Y+ n_gas] = (dYdz + (Y * sum_continuity
                                             - phi * np.multiply(wdot_g,W_g)
                                             - As * np.multiply(sdot_g,W_g)) / rhou)
    res[offset_Y + membrane_sp_ind] -= area2vol * membrane_flux * membrane_present

    # Conservation of site fractions (algebraic constraints in this example)
    res[offset_Z : offset_Z + n_surf] = sdot_s

    # For the species with largest site coverage (k_large), solve the constraint
    # equation: sum(Zk) = 1.
    # The residual function for 'k_large' would be 'res[k_large] = 1 - sum(Zk)'
    # Note that here sum(Zk) will be the sum of coverages for all surface species,
    # including the 'k_large' species.
    ind_large = np.argmax(coverages)
    res[offset_Z + ind_large] = 1 - np.sum(coverages)

    # Conservation of momentum
    u = rhou / rho
    res[offset_p] = dpdz + phi * mu * u / B_g

    # Conservation of energy
    res[offset_T] = dTdz - 0  # isothermal condition
    # Note: One can just not solve the energy equation by keeping temperature constant.
    # But since 'T' is used as the dependent variable, the residual is res[T] = dTdz - 0
    # One can also write res[T] = 0 directly, but that leads to a solver failure due to
    # singular jacobian

    if solve_energy:
        conv_term = (4 / D_h) * h_coeff * (T_wall - T) * (2 * np.pi * R)
        chem_term = np.sum(hk_g * (phi * wdot_g + As * sdot_g))
        res[offset_T] -= (conv_term - chem_term) / (rhou * cp)

# %%
# Calculate the spatial derivatives at the inlet that will be used as the initial
# conditions for the IDA solver

# Initialize yPrime to 0 and call residual to get initial derivatives
y0 = np.hstack((rhou0, p_in, T_in, Yk_0, Zk_0))
yprime0 = np.zeros(n_var)
res = np.zeros(n_var)
residual(0, y0, yprime0, res)
yprime0 = -res

# %%
# Solve the system of DAEs using IDA solver
# -----------------------------------------

solver = dae(
    "ida",
    residual,
    first_step_size=1e-15,
    atol=1e-14,  # absolute tolerance for solution
    rtol=1e-06,  # relative tolerance for solution
    algebraic_vars_idx=[np.arange(offset_Y + n_gas, offset_Z + n_surf, 1)],
    max_steps=8000,
    one_step_compute=True,
    old_api=False,  # forces use of new api (namedtuple)
)

distance = []
solution = []
state = solver.init_step(0.0, y0, yprime0)

# Note that here the variable t is an internal variable used in scikits. In this
# example, it represents the natural variable z, which corresponds to the axial distance
# inside the reactor.
while state.values.t < L:
    distance.append(state.values.t)
    solution.append(state.values.y)
    state = solver.step(L)

distance = np.array(distance)
solution = np.array(solution)
print(state)

# %%
# Plot results
# ------------

plt.rcParams['figure.constrained_layout.use'] = True

# %%
# Pressure and temperature
# ^^^^^^^^^^^^^^^^^^^^^^^^

# Plot gas pressure profile along the flow direction
f, ax = plt.subplots(figsize=(4,3))
ax.plot(distance, solution[:, offset_p], color="C0")
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Pressure (Pa)")

# Plot gas temperature profile along the flow direction
f, ax = plt.subplots(figsize=(4,3))
ax.plot(distance, solution[:, offset_T], color="C1")
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Temperature (K)")

# %%
# Mass fractions of the gas-phase species
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Plot major and minor gas species separately
minor_idx = []
major_idx = []
for j, name in enumerate(gas.species_names):
    mean = np.mean(solution[:, offset_Y + j])
    if mean <= 0.1:
        minor_idx.append(j)
    else:
        major_idx.append(j)

# Major gas-phase species
f, ax = plt.subplots(figsize=(4,3))
for j in major_idx:
    ax.plot(distance, solution[:, offset_Y + j], label=gas.species_names[j])
ax.legend(fontsize=12, loc="best")
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Mass Fraction")

# Minor gas-phase species
f, ax = plt.subplots(figsize=(4,3))
for j in minor_idx:
    ax.plot(distance, solution[:, offset_Y + j], label=gas.species_names[j])
ax.legend(fontsize=12, loc="best")
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Mass Fraction")

# %%
# Site fractions of the surface-phase species
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Plot major and minor surface species separately
minor_idx = []
major_idx = []
for j, name in enumerate(surf.species_names):
    mean = np.mean(solution[:, offset_Z + j])
    if mean <= 0.1:
        minor_idx.append(j)
    else:
        major_idx.append(j)

# Major surf-phase species
f, ax = plt.subplots(figsize=(4,3))
for j in major_idx:
    ax.plot(distance, solution[:, offset_Z + j], label=surf.species_names[j])
ax.legend(fontsize=12, loc="best")
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Site Fraction")

# Minor surf-phase species
f, ax = plt.subplots(figsize=(4,3))
for j in minor_idx:
    ax.plot(distance, solution[:, offset_Z + j], label=surf.species_names[j])
ax.legend(fontsize=12, loc="best")
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Site Fraction")

#%%
# References
# ----------
#
# .. [1] G. Kogekar (2021). "Computationally efficient and robust models of non-ideal
#        thermodynamics, gas-phase kinetics and heterogeneous catalysis in chemical
#        reactors", Doctoral dissertation.
#
# .. [2] B. Kee, C. Karakaya, H. Zhu, S. DeCaluwe, and R.J. Kee (2017). "The Influence
#        of Hydrogen-Permeable Membranes and Pressure on Methane Dehydroaromatization in
#        Packed-Bed Catalytic Reactors", *Industrial & Engineering Chemistry Research*
#        56, 13:3551-3559.
#
# .. [3] R.J. Kee, M.E. Coltrin, P. Glarborg, and H. Zhu (2018). *Chemically Reacting
#        Flow: Theory, Modeling and Simulation*, Wiley.
#
# .. [4] Z. Zhang, C. Karakaya, R.J. Kee, J. Douglas Way, C. Wolden (2019).
#        "Barium-Promoted Ruthenium Catalysts on Yittria-Stabilized Zirconia Supports
#        for Ammonia Synthesis", *ACS Sustainable Chemistry & Engineering*
#        7:18038-18047.

# sphinx_gallery_thumbnail_number = -1

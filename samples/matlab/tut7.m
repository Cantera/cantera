%  Tutorial 7: Thermodynamic Properties
%
help tut7

% A variety of thermodynamic property methods are provided.
gas = Air
set(gas,'T',800,'P',oneatm)

% temperature, pressure, density
T = temperature(gas)
P = pressure(gas)
rho = density(gas)
n = molarDensity(gas)

% species non-dimensional properties
hrt = enthalpies_RT(gas)            % vector of h_k/RT

% mixture properties per mole
hmole = enthalpy_mole(gas)
umole = intEnergy_mole(gas)
smole = entropy_mole(gas)
gmole = gibbs_mole(gas)

% mixture properties per unit mass
hmass = enthalpy_mass(gas)
umass = intEnergy_mass(gas)
smass = entropy_mass(gas)
gmass = gibbs_mass(gas)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
cleanup

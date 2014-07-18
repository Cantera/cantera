function dydt = conuv(t, y, gas, mw) %#ok<INUSL>
% CONUV ODE system for a constant-volume, adiabatic reactor.
%
%    Function CONUV evaluates the system of ordinary differential
%    equations for an adiabatic, constant-volume,
%    zero-dimensional reactor. It assumes that the 'gas' object
%    represents a reacting ideal gas mixture.


% Set the state of the gas, based on the current solution vector.
setMassFractions(gas, y(2:end), 'nonorm');
set(gas, 'T', y(1), 'Rho', density(gas));
nsp = nSpecies(gas);

% energy equation
wdot = netProdRates(gas);
tdot = - temperature(gas) * gasconstant * (enthalpies_RT(gas) - ones(nsp,1))' ...
    * wdot / (density(gas)*cv_mass(gas));

% set up column vector for dydt
dydt = [ tdot
    zeros(nsp, 1) ];

% species equations
rrho = 1.0/density(gas);
for i = 1:nsp
    dydt(i+1) = rrho*mw(i)*wdot(i);
end

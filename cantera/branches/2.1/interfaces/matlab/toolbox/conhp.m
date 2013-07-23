function dydt = conhp(t,y,gas,mw) %#ok<INUSL>
% CONHP ODE system for a constant-pressure, adiabatic reactor.
%
%    Function CONHP evaluates the system of ordinary differential
%    equations for an adiabatic, constant-pressure,
%    zero-dimensional reactor. It assumes that the 'gas' object
%    represents a reacting ideal gas mixture.


% Set the state of the gas, based on the current solution vector.
set(gas, 'T', y(1), 'P', pressure(gas), 'Y', y(2:end));

% energy equation
wdot = netProdRates(gas);
tdot = - temperature(gas) * gasconstant * enthalpies_RT(gas)' ...
    * wdot / (density(gas)*cp_mass(gas));

% set up column vector for dydt
dydt = [ tdot
    zeros(53,1) ];

% species equations
rrho = 1.0/density(gas);
nsp = nSpecies(gas);
for i = 1:nsp
    dydt(i+1) = rrho*mw(i)*wdot(i);
end

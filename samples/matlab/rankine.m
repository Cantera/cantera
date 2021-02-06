function [work, efficiency] = rankine(t1, p2, eta_pump, eta_turbine)
%
% This example computes the efficiency of a simple vapor power cycle.
%
help rankine

% create an object representing water
w = Water;

% start with saturated liquid water at t1
set(w,'T',t1,'Liquid',1.0);
p1 = pressure(w);

% pump it to p2
pump_work = pump(w, p2, eta_pump);
h2 = enthalpy_mass(w);
p2 = pressure(w);

% heat to saturated vapor
set(w,'P',p2,'Vapor',1.0);
h3 = enthalpy_mass(w);

heat_added = h3 - h2;

% expand adiabatically back to the initial pressure
work = expand(w, p1, eta_turbine);

% compute the efficiency
efficiency = (work - pump_work)/heat_added;


function w = pump(fluid, pfinal, eta)
% PUMP - Adiabatically pump a fluid to pressure pfinal, using a pump
% with isentropic efficiency eta.
%
h0 = enthalpy_mass(fluid);
s0 = entropy_mass(fluid);
set(fluid, 'S', s0, 'P', pfinal);
h1s = enthalpy_mass(fluid);
isentropic_work = h1s - h0;
actual_work = isentropic_work / eta;
h1 = h0 + actual_work;
set(fluid, 'H',h1, 'P',pfinal);
w = actual_work;


function w = expand(fluid, pfinal, eta)
% EXPAND - Adiabatically expand a fluid to pressure pfinal, using a
% turbine with isentropic efficiency eta.
%
h0 = enthalpy_mass(fluid);
s0 = entropy_mass(fluid);
set(fluid, 'S', s0, 'P', pfinal);
h1s = enthalpy_mass(fluid);
isentropic_work = h0 - h1s;
actual_work = isentropic_work * eta;
h1 = h0 - actual_work;
set(fluid, 'H',h1, 'P',pfinal);
w = actual_work;

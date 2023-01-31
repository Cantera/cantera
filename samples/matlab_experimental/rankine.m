% RANKINE - This example computes the efficiency of a simple vapor power cycle.
%
% Keywords: thermodynamics, thermodynamic cycle, non-ideal fluid

clear all
close all

help rankine

% Initialize parameters
eta_pump = 0.6;
eta_turbine = 0.8;
p_max = 8.0 * OneAtm;
t1 = 300.0;

% create an object representing water
w = Water;

% start with saturated liquid water at t1
basis = 'mass';
w.setState_Tsat(t1, 1.0);
h1 = w.H;
p1 = w.P;

% pump it to p2
pump_work = pump(w, p_max, eta_pump);
h2 = w.H;

% heat to saturated vapor
w.setState_Psat(p_max, 1.0);
h3 = w.H;

heat_added = h3 - h2;

% expand adiabatically back to the initial pressure
turbine_work = expand(w, p1, eta_turbine);

% compute the efficiency
efficiency = (turbine_work - pump_work) / heat_added;
disp(sprintf('efficiency = %d', efficiency));

function work = pump(fluid, pfinal, eta)
    % PUMP - Adiabatically pump a fluid to pressure pfinal, using a pump
    % with isentropic efficiency eta.

    fluid.basis = 'mass';
    h0 = fluid.H;
    s0 = fluid.S;
    fluid.SP = {s0, pfinal};
    h1s = fluid.H;
    isentropic_work = h1s - h0;
    actual_work = isentropic_work / eta;
    h1 = h0 + actual_work;
    fluid.HP = {h1, pfinal};
    work = actual_work;
end

function work = expand(fluid, pfinal, eta)
    % EXPAND - Adiabatically expand a fluid to pressure pfinal, using a
    % turbine with isentropic efficiency eta.

    fluid.basis = 'mass';
    h0 = fluid.H;
    s0 = fluid.S;
    fluid.SP = {s0, pfinal};
    h1s = fluid.H;
    isentropic_work = h0 - h1s;
    actual_work = isentropic_work * eta;
    h1 = h0 - actual_work;
    fluid.HP = {h1, pfinal};
    work = actual_work;
end

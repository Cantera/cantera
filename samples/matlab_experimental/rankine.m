function [work, efficiency] = rankine(t1, p2, eta_pump, eta_turbine)
    % This example computes the efficiency of a simple vapor power cycle.
    %
    % Keywords: thermodynamics, thermodynamic cycle, non-ideal fluid

    help rankine

    % create an object representing water
    w = Water;

    % start with saturated liquid water at t1
    w.setState_Tsat(t1, 1.0);
    p1 = w.P;

    % pump it to p2
    basis = 'mass';
    pump_work = pump(w, p2, eta_pump);
    h2 = w.H;
    p2 = w.P;

    % heat to saturated vapor
    w.setState_Psat(p2, 1.0);
    h3 = w.H;

    heat_added = h3 - h2;

    % expand adiabatically back to the initial pressure
    work = expand(w, p1, eta_turbine);

    % compute the efficiency
    efficiency = (work - pump_work)/heat_added;
end


function w = pump(fluid, pfinal, eta)
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
    w = actual_work;
end


function w = expand(fluid, pfinal, eta)
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
    w = actual_work;
end

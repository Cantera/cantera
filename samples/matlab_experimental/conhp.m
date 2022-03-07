function dydt = conhp(t, y, gas, mw) %#ok<INUSL>
    % CONHP - ODE system for a constant-pressure, adiabatic reactor.
    %
    %    Function CONHP evaluates the system of ordinary differential
    %    equations for an adiabatic, constant-pressure,
    %    zero-dimensional reactor. It assumes that the 'gas' object
    %    represents a reacting ideal gas mixture.

    % Set the state of the gas, based on the current solution vector.
    gas.Y = y(2: end);
    gas.TP = {y(1), gas.P};
    nsp = gas.nSpecies;

    % energy equation
    wdot = gas.netProdRates;
    H = gas.enthalpies_RT';
    gas.basis = 'mass';
    tdot = - gas.T * gasconstant /(gas.D * gas.cp).* wdot * H;

    % set up column vector for dydt
    dydt = [tdot
            zeros(nsp, 1)];

    % species equations
    rrho = 1.0/gas.D;
    for i = 1:nsp
        dydt(i+1) = rrho * mw(i) * wdot(i);
    end
end

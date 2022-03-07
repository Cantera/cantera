function dydt = conuv(t, y, gas, mw) %#ok<INUSL>
    % CONUV ODE system for a constant-volume, adiabatic reactor.
    %
    %    Function CONUV evaluates the system of ordinary differential
    %    equations for an adiabatic, constant-volume,
    %    zero-dimensional reactor. It assumes that the 'gas' object
    %    represents a reacting ideal gas mixture.


    % Set the state of the gas, based on the current solution vector.
    gas.Y = y(2:end);
    gas.TD = {y(1), gas.D};
    nsp = gas.nSpecies;

    % energy equation
    wdot = gas.netProdRates;
    H = (gas.enthalpies_RT - ones(1, nsp))';
    gas.basis = 'mass';
    tdot = - gas.T * gasconstant / (gas.D * gas.cv).* wdot * H;

    % set up column vector for dydt
    dydt = [tdot
            zeros(nsp, 1)];

    % species equations
    rrho = 1.0/gas.D;
    for i = 1:nsp
        dydt(i+1) = rrho * mw(i) * wdot(i);
    end

end

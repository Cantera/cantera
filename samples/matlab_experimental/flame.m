function f = flame(gas, left, flow, right)
    %% Utility for flame setup
    %
    % Used by the :doc:`flame1.m <flame1>` and :doc:`flame2.m <flame2>` examples.

    % Check input parameters
    if nargin ~= 4
        error('wrong number of input arguments.');
    end

    if ~gas.isIdealGas
        error('gas object must represent an ideal gas mixture.');
    end

    if ~isa(left, 'Inlet')
        error('burner object of wrong type.');
    end

    if ~isa(flow, 'Flow1D')
        error('flow object of wrong type.');
    end

    flametype = 0;

    if isa(right, 'ReactingSurface')
        flametype = 1;
    elseif isa(right, 'Inlet')
        flametype = 3;
    end

    % create the container object
    f = Sim1D({left flow right});

    rho0 = gas.massDensity;

    % find the adiabatic flame temperature and corresponding
    % equilibrium composition
    gas.equilibrate('HP');
    teq = gas.T;
    yeq = gas.Y;

    z1 = 0.2;
    mdot0 = left.massFlux;
    mdot1 = right.massFlux;
    t0 = left.T;
    u0 = mdot0 / rho0;

    if flametype == 0
        t1 = teq;
        u1 = -mdot1 / gas.massDensity;
    else
        t1 = right.T;
        u1 = -mdot1 / rho0;
    end

    flow.setProfile('velocity', [0.0, 1.0], [u0, u1]);
    flow.setProfile('spreadRate', [0.0, 1.0], [0.0, 0.0]);
    flow.setProfile('T', [0.0, z1, 1.0], [t0, teq, t1]);

    for n = 1:gas.nSpecies
        nm = gas.speciesName(n);

        if flametype == 3
            y1 = right.massFraction(n);
        else
            y1 = yeq(n);
        end

        flow.setProfile(nm{:}, [0.0, z1, 1.0], [left.massFraction(n), yeq(n), y1]);
    end

    % set minimal grid refinement criteria
    flow.setRefineCriteria(10.0, 0.8, 0.8);
end

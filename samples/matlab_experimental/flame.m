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

    % set default initial profiles.
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

    if flametype == 0
        t1 = teq;
        mdot1 = -mdot0;
    else
        t1 = right.T;
    end

    f.setProfile(2, {'velocity', 'spread_rate'}, [0.0 1.0
                                                  mdot0 / rho0 -mdot1 / rho0
                                                  0.0 0.0]);
    f.setProfile(2, 'T', [0.0, 1.0
                          t0, t1]);

    for n = 1:gas.nSpecies
        nm = gas.speciesName(n);

        if strcmp(nm, 'H') || strcmp(nm, 'OH') || strcmp(nm, 'O') || ...
            strcmp(nm, 'HO2')
            yint = 1.0 * yeq(n);
        else
            yint = yeq(n);
        end

        if flametype == 3
            y1 = right.massFraction(n);
        else
            y1 = yeq(n);
        end

        f.setProfile(2, nm, [0, 1.0
                    left.massFraction(n), y1]);
    end

    % set minimal grid refinement criteria
    f.setRefineCriteria(2, 10.0, 0.8, 0.8);
end

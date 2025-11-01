function ignite_hp(gas)
    %% Constant pressure ignition with user-specified equations
    %
    % Solves the same ignition problem as :doc:`reactor1.m <reactor1>`, but uses
    % the local function ``REACTOR_ODE`` to implement the governing equations for
    % an adiabatic, constant-pressure, zero-dimensional reactor.
    %
    % Requires: cantera >= 3.2.0
    %
    % .. tags:: Matlab, combustion, user-defined model, ignition delay, plotting

    tic
    help ignite_hp

    if nargin == 0
        gas = Solution('gri30.yaml', 'gri30');
    end

    mw = gas.molecularWeights;
    gas.TPX = {1001.0, OneAtm, 'H2:2,O2:1,N2:4'};
    gas.basis = 'mass';
    nsp = gas.nSpecies;

    y0 = [gas.T
          gas.Y'];
    tel = [0, 0.001];
    options = odeset('RelTol', 1.e-5, 'AbsTol', 1.e-12, 'Stats', 'on');
    t0 = cputime;
    out = ode15s(@reactor_ode, tel, y0, options, gas, mw, nsp);
    disp(['CPU time = ' num2str(cputime - t0)]);

    if nargout == 0
        % plot the temperature and OH mole fractions.
        figure(1);
        plot(out.x, out.y(1, :));
        xlabel('time');
        ylabel('Temperature');
        title(['Final T = ' num2str(out.y(1, end)), ' K']);

        figure(2);
        ioh = gas.speciesIndex('OH');
        plot(out.x, out.y(1 + ioh, :));
        xlabel('time');
        ylabel('Mass Fraction');
        title('OH Mass Fraction');
    end

    toc
end

%%
% Local Functions
% ---------------

function dydt = reactor_ode(t, y, gas, mw, nsp)
    % ODE system for a constant-pressure, adiabatic reactor
    %
    % Function ``REACTOR_ODE`` evaluates the system of ordinary differential equations
    % for an adiabatic, constant-pressure, zero-dimensional reactor.
    % It assumes that the ``gas`` object represents a reacting ideal gas mixture.

    % Set the state of the gas, based on the current solution vector.
    gas.TPY = {y(1), gas.P, y(2:end)};

    % energy equation
    wdot = gas.netProdRates;
    H = gas.partialMolarEnthalpies';
    tdot =- 1 / (gas.D * gas.cp) .* wdot * H;

    % set up column vector for dydt
    dydt = [tdot
            zeros(nsp, 1)];

    % species equations
    rrho = 1.0 / gas.D;

    for i = 1:nsp
        dydt(i + 1) = rrho * mw(i) * wdot(i);
    end

end

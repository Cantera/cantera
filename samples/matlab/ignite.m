function ignite(g)
    %% Adiabatic, constant pressure reactor
    %
    % This example solves the same problem as :doc:`reactor1.m <reactor1>`, but does it
    % using one of MATLAB's ODE integrators, rather than using the Cantera Reactor
    % class. The governing equations are implemented using the local function
    % ``REACTOR_ODE``.
    %
    % Requires: cantera >= 3.2.0
    %
    % .. tags:: Matlab, combustion, reactor network, ignition delay, plotting

    tic
    help ignite

    if nargin == 1
        gas = g;
    else
        gas = Solution('gri30.yaml', 'gri30');
    end

    % set the initial conditions

    gas.TPX = {1001.0, OneAtm, 'H2:2,O2:1,N2:4'};
    gas.basis = 'mass';
    y0 = [gas.U, 1.0 / gas.massDensity, gas.Y];

    time_interval = [0, 0.001];
    options = odeset('RelTol', 1.e-5, 'AbsTol', 1.e-12, 'Stats', 'on');

    t0 = cputime;
    out = ode15s(@reactor_ode, time_interval, y0, options, gas, ...
                 @vdot, @area, @heatflux);
    disp(['CPU time = ' num2str(cputime - t0)]);

    toc
end

%%
% Local Functions
% ---------------

function dydt = reactor_ode(t, y, gas, vdot, area, heatflux)
    % ODE system for a generic zero-dimensional reactor
    %
    % Function ``REACTOR_ODE`` evaluates the system of ordinary differential equations
    % for a zero-dimensional reactor with arbitrary heat transfer and volume change.
    % Used in :doc:`ignite.m <ignite>`.
    %
    % Solution vector components:
    %    :y(1):     Total internal energy U
    %    :y(2):     Volume V
    %    :y(3):     Mass of species 1
    %    :....:
    %    :y(nsp+2): Mass of last species
    %

    [m, n] = size(y);
    dydt = zeros(m, n);

    for j = 1:n
        this_y = y(:, j);
        int_energy = this_y(1);
        vol = this_y(2);
        masses = this_y(3:end);

        % evaluate the total mass, and the specific internal energy and volume.
        total_mass = sum(masses);
        u_mass = int_energy / total_mass;
        v_mass = vol / total_mass;

        % set the state of the gas by specifying (u,v,{Y_k})
        gas.UVY = {u_mass, v_mass, masses};
        p = gas.P;

        % volume equation
        vdt = feval(vdot, t, vol, gas);

        % energy equation
        a = feval(area, t, vol);
        q = feval(heatflux, t, gas);
        udt = -p * vdt + a * q;

        % species equations
        k = gas.netProdRates;
        rho_inv = 1 / gas.massDensity;
        MW = gas.molecularWeights;
        massProdRate = rho_inv .* k .* MW;
        ydt = total_mass * massProdRate;

        % set up column vector for dydt
        dydt(:, j) = [udt, vdt, ydt];
    end

end


function v = vdot(t, vol, gas)
    % Time-varying boundary conditions.
    %
    % The functions below may be defined arbitrarily to set the reactor
    % boundary conditions - the rate of change of volume, the heat
    % flux, and the area.
    %
    % Rate of change of volume. Any arbitrary function may be implemented.
    %
    % Input arguments:
    %    :t:      time
    %    :vol:    volume
    %    :gas:    ideal gas object

    %v = 0.0;                                 %uncomment for constant volume
    v = 1.e11 * (gas.P - 101325.0); % holds pressure very
    % close to 1 atm
end


function q = heatflux(t, gas)
    % heat flux (W/m^2).
    q = 0.0; % adiabatic
end


function a = area(t, vol)
    % surface area (m^2). Used only to compute heat transfer.
    a = 1.0;
end


function pv = output(s, gas)
    % Since the solution variables used by the ``reactor`` function are
    % not necessarily those desired for output, this function is called
    % after the integration is complete to generate the desired
    % outputs.

    times = s.x;
    soln = s.y;
    [~, n] = size(times);
    pv = zeros(gas.nSpecies + 4, n);

    gas.TP = {1001.0, OneAtm};

    for j = 1:n
        ss = soln(:, j);
        y = ss(3:end);
        mass = sum(y);
        u_mass = ss(1) / mass;
        v_mass = ss(2) / mass;
        gas.UVY = {u_mass, v_mass, y};

        pv(1, j) = times(j);
        pv(2, j) = gas.T;
        pv(3, j) = gas.D;
        pv(4, j) = gas.P;
        pv(5:end, j) = y;
    end

    % plot the temperature and OH mass fractions.
    clf;

    subplot(1, 2, 1);
    plot(pv(1, :), pv(2, :));
    xlabel('time');
    ylabel('Temperature');
    title(['Final T = ' num2str(pv(2, end)) ' K']);

    subplot(1, 2, 2)
    ioh = gas.speciesIndex('OH');
    plot(pv(1, :), pv(4 + ioh, :));
    xlabel('time');
    ylabel('Mass Fraction');
    title('OH Mass Fraction');
end

function plotdata = ignite(g)
    %% Adiabatic, constant pressure reactor
    %
    % This example solves the same problem as :doc:`reactor1.m <reactor1>`, but does it
    % using one of MATLAB's ODE integrators, rather than using the Cantera Reactor
    % class. See :doc:`reactor_ode.m <reactor_ode>` for the implementation of the
    % governing equations.
    %
    % .. tags:: Matlab, combustion, reactor network, ignition delay, plotting

    clear all
    close all

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
    y0 = [gas.U
          1.0 / gas.D
          gas.Y'];

    time_interval = [0 0.001];
    options = odeset('RelTol', 1.e-5, 'AbsTol', 1.e-12, 'Stats', 'on');

    t0 = cputime;
    out = ode15s(@reactor_ode, time_interval, y0, options, gas, ...
                @vdot, @area, @heatflux);
    disp(['CPU time = ' num2str(cputime - t0)]);

    plotdata = output(out, gas);

    %% Time-varying boundary conditions
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
    function v = vdot(t, vol, gas)
        %v = 0.0;                                 %uncomment for constant volume
        v = 1.e11 * (gas.P - 101325.0); % holds pressure very
        % close to 1 atm
    end

    %%
    % heat flux (W/m^2).
    function q = heatflux(t, gas)
        q = 0.0; % adiabatic
    end

    %%
    % surface area (m^2). Used only to compute heat transfer.
    function a = area(t, vol)
        a = 1.0;
    end

    %%
    % Since the solution variables used by the ``reactor`` function are
    % not necessarily those desired for output, this function is called
    % after the integration is complete to generate the desired
    % outputs.

    function pv = output(s, gas)
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
            gas.Y = y;
            gas.UV = {u_mass, v_mass};

            pv(1, j) = times(j);
            pv(2, j) = gas.T;
            pv(3, j) = gas.D;
            pv(4, j) = gas.P;
            pv(5:end, j) = y;
        end

        % plot the temperature and OH mass fractions.
        figure(1);
        plot(pv(1, :), pv(2, :));
        xlabel('time');
        ylabel('Temperature');
        title(['Final T = ' num2str(pv(2, end)) ' K']);

        figure(2);
        ioh = gas.speciesIndex('OH');
        plot(pv(1, :), pv(4 + ioh, :));
        xlabel('time');
        ylabel('Mass Fraction');
        title('OH Mass Fraction');
    end

    toc
end

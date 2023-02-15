function reactor1(g)
    %% REACTOR1 Zero-dimensional kinetics: adiabatc, constant pressure.
    %
    % This example illustrates how to use class 'Reactor' for zero-dimensional
    % kinetics simulations. Here the parameters are set so that the reactor is
    % adiabatic and very close to constant pressure.
    %
    % Keywords: combustion, reactor network, ignition delay, plotting

    clear all
    close all

    tic
    help reactor1

    if nargin == 1
        gas = g;
    else
        gas = Solution('gri30.yaml', 'gri30', 'None');
    end

    P = 101325.0;
    % set the initial conditions
    gas.TP = {1001.0, P};
    nsp = gas.nSpecies;
    xx = zeros(nsp, 1);
    xx(1) = 0.285;
    xx(4) = 0.142;
    xx(48) = 0.573;
    gas.X = xx;

    % create a reactor, and insert the gas
    r = IdealGasReactor(gas);

    % create a reservoir to represent the environment
    a = Solution('air.yaml', 'air', 'None');
    a.TP = {a.T, P};
    env = Reservoir(a);

    % Define a wall between the reactor and the environment and
    % make it flexible, so that the pressure in the reactor is held
    % at the environment pressure.
    w = Wall(r, env);

    % set expansion parameter. dV/dt = KA(P_1 - P_2)
    w.expansionRateCoeff = 1.0e6;

    % set wall area
    w.area = 1.0;

    % create a reactor network and insert the reactor:
    network = ReactorNet({r});

    nSteps = 100;
    tim(nSteps) = 0;
    temp(nSteps) = 0;
    x(nSteps, 3) = 0;
    t = 0.0;
    dt = 1.0e-5;
    t0 = cputime;

    for n = 1:nSteps
        t = t + dt;
        network.advance(t);
        tim(n) = network.time;
        temp(n) = r.T;
        x(n, 1:3) = gas.moleFraction({'OH', 'H', 'H2'});
    end

    disp(['CPU time = ' num2str(cputime - t0)]);

    clf;
    subplot(2, 2, 1);
    plot(tim, temp);
    xlabel('Time (s)');
    ylabel('Temperature (K)');
    subplot(2, 2, 2)
    plot(tim, x(:, 1));
    xlabel('Time (s)');
    ylabel('OH Mole Fraction (K)');
    subplot(2, 2, 3)
    plot(tim, x(:, 2));
    xlabel('Time (s)');
    ylabel('H Mole Fraction (K)');
    subplot(2, 2, 4)
    plot(tim, x(:, 3));
    xlabel('Time (s)');
    ylabel('H2 Mole Fraction (K)');

    toc
end

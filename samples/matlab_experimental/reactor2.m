function reactor2(g)
    %% REACTOR2 - Zero-dimensional kinetics: adiabatic, constant volume.
    %
    % This example illustrates how to use class 'Reactor' for zero-dimensional
    % kinetics simulations. Here the parameters are set so that the reactor is
    % adiabatic and constant volume.
    %
    % Keywords: combustion, reactor network, ignition delay, plotting

    clear all
    close all

    tic
    help reactor2

    if nargin == 1
        gas = g;
    else
        gas = Solution('gri30.yaml', 'gri30', 'None');
    end

    % set the initial conditions
    gas.TPX = {1001.0, oneatm, 'H2:2,O2:1,N2:4'};

    % create a reactor, and insert the gas
    r = IdealGasReactor(gas);

    % create a reactor network and insert the reactor
    network = ReactorNet({r});

    nSteps = 100;
    tim(nSteps) = 0;
    temp(nSteps) = 0;
    x(nSteps, 3) = 0;
    t = 0;
    dt = 1.0e-5;
    t0 = cputime;

    for n = 1:100
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

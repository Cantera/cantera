function periodic_cstr
%
%  A CSTR with steady inputs but periodic interior state.
%
%  A stoichiometric hydrogen/oxygen mixture is introduced and reacts to produce
%  water. But since water has a large efficiency as a third body in the chain
%  termination reaction
%
%         H + O2 + M = HO2 + M
%
%  as soon as a significant amount of water is produced the reaction
%  stops. After enough time has passed that the water is exhausted from
%  the reactor, the mixture explodes again and the process
%  repeats. This explanation can be verified by decreasing the rate for
%  reaction 7 in file 'h2o2.yaml' and re-running the example.
%
%  Acknowledgments: The idea for this example and an estimate of the
%  conditions needed to see the oscillations came from Bob Kee,
%  Colorado School of Mines
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    help periodic_cstr

    % create the gas mixture
    gas = Solution('h2o2.yaml');

    % pressure = 60 Torr, T = 770 K
    p = 60.0 * 133.3;
    t = 770.0;

    gas.TPX = {300, p, 'H2:2, O2:1'};

    % create an upstream reservoir that will supply the reactor.  The
    % temperature, pressure, and composition of the upstream reservoir are
    % set to those of the 'gas' object at the time the reservoir is
    % created.
    upstream = Reservoir(gas);

    % Now set the gas to the initial temperature of the reactor, and create
    % the reactor object.
    gas.TP = {t, p};
    cstr = IdealGasReactor(gas);

    % Set its volume to 10 cm^3. In this problem, the reactor volume is
    % fixed, so the initial volume is the volume at all later times.
    cstr.setInitialVolume(10.0*1.0e-6);

    % We need to have heat loss to see the oscillations. Create a
    % reservoir to represent the environment, and initialize its
    % temperature to the reactor temperature.
    env = Reservoir(gas);

    % Create a heat-conducting wall between the reactor and the
    % environment. Set its area, and its overall heat transfer
    % coefficient. Larger U causes the reactor to be closer to isothermal.
    % If U is too small, the gas ignites, and the temperature spikes and
    % stays high.
    w = Wall;
    w.install(cstr, env);
    w.area = 1.0;
    w.setHeatTransferCoeff(0.02);

    % Connect the upstream reservoir to the reactor with a mass flow
    % controller (constant mdot). Set the mass flow rate to 1.25 sccm.
    sccm = 1.25;
    vdot = sccm * 1.0e-6/60.0 * ((oneatm / gas.P) * ( gas.T / 273.15));  % m^3/s
    mdot = gas.D * vdot;   % kg/s
    mfc = MassFlowController;
    mfc.install(upstream, cstr);
    mfc.setMassFlowRate(mdot);

    % now create a downstream reservoir to exhaust into.
    downstream = Reservoir(gas);

    % connect the reactor to the downstream reservoir with a valve, and
    % set the coefficient sufficiently large to keep the reactor pressure
    % close to the downstream pressure of 60 Torr.
    v = Valve;
    v.install(cstr, downstream);
    v.setValveCoeff(1.0e-9);

    % create the network
    network = ReactorNet({cstr});

    % now integrate in time
    tme = 0.0;
    dt   = 0.1;

    n = 0;
    while tme < 300.0
        n = n + 1;
        tme = tme + dt;
        network.advance(tme);
        tm(n) = tme;
        y(1,n) = cstr.massFraction('H2');
        y(2,n) = cstr.massFraction('O2');
        y(3,n) = cstr.massFraction('H2O');
    end
    clf
    figure(1)
    plot(tm,y)
    legend('H2','O2','H2O')
    title('Mass Fractions')
end

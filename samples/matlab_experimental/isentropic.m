function isentropic(g)
    %% Isentropic, adiabatic flow
    %
    % In this example, the area ratio vs. Mach number curve is computed for a
    % hydrogen/nitrogen gas mixture.
    %
    % .. tags:: Matlab, thermodynamics, compressible flow, plotting

    clear all
    close all

    tic
    help isentropic

    if nargin == 1
        gas = g;
    else
        gas = Solution('gri30.yaml', 'gri30');
    end

    % set the stagnation state
    gas.TPX = {1200.0, 10.0 * OneAtm, 'H2:1,N2:0.1'};
    gas.basis = 'mass';
    s0 = gas.S;
    h0 = gas.H;
    p0 = gas.P;

    mdot = 1; % arbitrary

    mach = [];
    a = [];
    i = 1;
    amin = 1.e14;

    % compute values for a range of pressure ratios
    for r = 0.005:0.0025:0.995
        p = p0 * r;

        % set the state using (p,s0)
        gas.SP = {s0, p};

        h = gas.H;
        rho = gas.D;

        v2 = 2.0 * (h0 - h); %   h + V^2/2 = h0
        v = sqrt(v2);
        a(i) = mdot / (rho * v); %   rho*v*A = constant

        if a(i) < amin
            amin = a(i);
        end

        mach(i) = v / gas.soundSpeed;
        i = i + 1;
    end

    a = a / amin;

    % plot results

    clf;
    plot(mach, a);
    ylabel('Area Ratio');
    xlabel('Mach Number');
    title('Isentropic Flow: Area Ratio vs. Mach Number');

    toc
end

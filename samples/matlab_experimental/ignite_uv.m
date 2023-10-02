function ignite_uv(gas)
    %% Constant volume ignition with user-specified equations
    %
    % Solves the same ignition problem as :doc:`reactor2.m <reactor2>`, except using
    % function :doc:`conuv.m <conuv>` to implement the governing equations.
    %
    % .. tags:: Matlab, combustion, user-defined model, ignition delay, plotting

    clear all
    close all

    tic
    help ignite_uv

    if nargin == 0
        gas = Solution('gri30.yaml', 'gri30');
    end

    mw = gas.molecularWeights;
    gas.TPX = {1001.0, OneAtm, 'H2:2,O2:1,N2:4'};

    y0 = [gas.T
          gas.X'];
    tel = [0, 0.001];
    options = odeset('RelTol', 1.e-5, 'AbsTol', 1.e-12, 'Stats', 'on');
    t0 = cputime;
    out = ode15s(@conuv, tel, y0, options, gas, mw);
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

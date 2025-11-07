%% Burner-stabilized flat flame
%
% This script simulates a burner-stabilized lean hydrogen-oxygen flame
% at low pressure. This example is equivalent to the Python
% :doc:`burner_flame.py <../python/onedim/burner_flame>` example.
%
% Requires: cantera >= 3.2.0
%
% .. tags:: Matlab, combustion, 1D flow, burner-stabilized flame, plotting

%%
% Problem Definition
% ------------------
%

tic
help burner_flame

t0 = cputime;  % record the starting time

%%
% **Set parameter values**

p = 0.05 * OneAtm;  % pressure
Tburner = 373.0;  % burner temperature
mdot = 0.06;  % kg/m^2/s

rxnmech = 'h2o2.yaml';  % reaction mechanism file
comp = 'H2:1.5, O2:1, AR:7';  % premixed gas composition

width = 0.5;
nz = 11;

logLevel = 1;  % amount of diagnostic output (0 to 5)
refineGrid = 1;  % 1 to enable refinement, 0 to disable
maxJacobianAge = [5, 10];

%%
% **Create the gas object**
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties

gas = Solution(rxnmech, 'ohmech', 'mixture-averaged');

% set its state to that of the unburned gas at the burner
gas.TPX = {Tburner, p, comp};
rhoIn = gas.massDensity;
Yin = gas.Y;

gas.equilibrate('HP');
rhoOut = gas.massDensity;
Yout = gas.Y;
Tad = gas.T;

%%
% **Create the flow object**

flow = UnstrainedFlow(gas, 'flow');
flow.P = p;
flow.setupUniformGrid(nz, width, 0.0);
flow.energyEnabled = false;

%%
% **Create the burner**
%
% The burner is an ``Inlet`` object. The temperature, mass flux,
% and composition (relative molar) may be specified.

burner = Inlet(gas, 'burner');
burner.T = Tburner;
burner.X = comp;
burner.massFlux = mdot;
uIn = mdot / rhoIn;

%%
% **Create the outlet**
%
% The type of flame is determined by the object that terminates
% the domain. An ``Outlet`` object imposes zero gradient boundary
% conditions for the temperature and mass fractions, and zero
% radial velocity and radial pressure gradient.

out = Outlet(gas, 'out');
uOut = mdot / rhoOut;

%%
% **Create the flame object**
%
% Once the component parts have been created, they can be assembled
% to create the flame object.
fl = Sim1D({burner, flow, out});
fl.setMaxJacAge(maxJacobianAge(1), maxJacobianAge(2));

% Supply initial guess
locs = [0.0, 0.2, 1.0];
flow.setProfile('velocity', locs, [uIn, uOut, uOut]);
flow.setProfile('T', locs, [Tburner, Tad, Tad]);

names = gas.speciesNames;
for i = 1:gas.nSpecies
    flow.setProfile(names{i}, locs, [Yin(i), Yout(i), Yout(i)]);
end

%%
% Solution
% --------
%
% Start with energy equation disabled

fl.solve(logLevel, refineGrid);

%%
% **Enable the energy equation**
%
% The energy equation will now be solved to compute the
% temperature profile. We also tighten the grid refinement
% criteria to get an accurate final solution.

flow.energyEnabled = true;
flow.setRefineCriteria(3.0, 0.05, 0.1);
fl.solve(logLevel, refineGrid);

%%
% **Show statistics and display results**
%
% Show statistics

fl.writeStats;
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);

toc

%%
% Plot results

clf;
subplot(2, 2, 1);
plotSolution(flow, 'T', 'Temperature [K]');
subplot(2, 2, 2);
plotSolution(flow, 'velocity', 'Axial Velocity [m/s]');
subplot(2, 2, 3);
plotSolution(flow, 'H2O', 'H2O Mass Fraction');
subplot(2, 2, 4);
plotSolution(flow, 'O2', 'O2 Mass Fraction');

%%
% Plotting Utility
% ----------------
%
% Local helper function to create plots

function plotSolution(domain, component, titleStr)
    % Utility for plotting a specific solution component
    z = domain.grid;
    x = domain.values(component);
    plot(z, x);
    xlabel('z (m)');
    ylabel(component);
    title(titleStr);
end

%% Burner-stabilized flat flame
%
% This script simulates a burner-stablized lean hydrogen-oxygen flame
% at low pressure.
%
% .. tags:: Matlab, combustion, 1D flow, burner-stabilized flame, plotting

%% Initialization

tic
help flame1

t0 = cputime; % record the starting time

%% Set parameter values

p = 0.05 * OneAtm; % pressure
tburner = 373.0; % burner temperature
mdot = 0.06; % kg/m^2/s

rxnmech = 'h2o2.yaml'; % reaction mechanism file
comp = 'H2:1.5, O2:1, AR:7'; % premixed gas composition

width = 0.5;
nz = 11;

logLevel = 1; % amount of diagnostic output (0 to 5)
refineGrid = 1; % 1 to enable refinement, 0 to disable
maxJacobianAge = [5, 10];

%% Create the gas object
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties

gas = Solution(rxnmech, 'ohmech', 'mixture-averaged');

% set its state to that of the unburned gas at the burner
gas.TPX = {tburner, p, comp};

%% Create the flow object

f = AxisymmetricFlow(gas, 'flow');
f.P = p;
f.setupUniformGrid(nz, width, 0.0);

%% Create the burner
%
%  The burner is an ``Inlet`` object. The temperature, mass flux,
%  and composition (relative molar) may be specified.

burner = Inlet(gas, 'burner');
burner.T = tburner;
burner.massFlux = mdot;
burner.X = comp;

%% Create the outlet
%
%  The type of flame is determined by the object that terminates
%  the domain. An ``Outlet`` object imposes zero gradient boundary
%  conditions for the temperature and mass fractions, and zero
%  radial velocity and radial pressure gradient.

s = Outlet(gas, 'out');

%% Create the flame object
%
% Once the component parts have been created, they can be assembled
% to create the flame object (see :doc:`flame.m <flame>`).
fl = flame(gas, burner, f, s);
fl.setMaxJacAge(maxJacobianAge(1), maxJacobianAge(2));

%%
% if the starting solution is to be read from a previously-saved
% solution, uncomment this line and edit the file name and solution id.

%restore(fl,'h2flame2.yaml', 'energy')
fl.solve(logLevel, refineGrid);

%% Enable the energy equation
%
%  The energy equation will now be solved to compute the
%  temperature profile. We also tighten the grid refinement
%  criteria to get an accurate final solution.

f.energyEnabled = true;
f.setRefineCriteria(3.0, 0.05, 0.1);
fl.solve(logLevel, refineGrid);

%% Show statistics

fl.writeStats;
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);

%% Make plots

clf;
subplot(2, 2, 1);
plotSolution(f, 'T');
title('Temperature [K]');
subplot(2, 2, 2);
plotSolution(f, 'velocity');
title('Axial Velocity [m/s]');
subplot(2, 2, 3);
plotSolution(f, 'H2O');
title('H2O Mass Fraction');
subplot(2, 2, 4);
plotSolution(f, 'O2');
title('O2 Mass Fraction');

toc

%% FLAME1 - A burner-stabilized flat flame
%
% This script simulates a burner-stablized lean hydrogen-oxygen flame
% at low pressure.
%
% Keywords: combustion, 1D flow, burner-stabilized flame, plotting

%% Initialization

clear all
close all
clc

tic
help flame1

t0 = cputime; % record the starting time

%% Set parameter values

p = 0.05 * oneatm; % pressure
tburner = 373.0; % burner temperature
mdot = 0.06; % kg/m^2/s

rxnmech = 'h2o2.yaml'; % reaction mechanism file
comp = 'H2:1.8, O2:1, AR:7'; % premixed gas composition

initial_grid = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.4, ...
                0.49, 0.5]; % m

tol_ss = {1.0e-5, 1.0e-13}; % {rtol atol} for steady-state problem
tol_ts = {1.0e-4, 1.0e-9}; % {rtol atol} for time stepping

loglevel = 1; % amount of diagnostic output (0
% to 5)

refine_grid = 1; % 1 to enable refinement, 0 to
% disable
max_jacobian_age = [5, 10];

%% Create the gas object
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties

gas = Solution(rxnmech, 'ohmech', 'Mix');

% set its state to that of the unburned gas at the burner
gas.TPX = {tburner, p, comp};

%% Create the flow object

f = AxisymmetricFlow(gas, 'flow');
f.P = p;
f.setupGrid(initial_grid);
f.setSteadyTolerances('default', tol_ss{:});
f.setTransientTolerances('default', tol_ts{:});

%% Create the burner
%
%  The burner is an Inlet object. The temperature, mass flux,
%  and composition (relative molar) may be specified.
burner = Inlet('burner');
burner.T = tburner;
burner.setMassFlowRate(mdot);
burner.setMoleFractions(comp);

%% Create the outlet
%
%  The type of flame is determined by the object that terminates
%  the domain. An Outlet object imposes zero gradient boundary
%  conditions for the temperature and mass fractions, and zero
%  radial velocity and radial pressure gradient.

s = Outlet('out');

%% Create the flame object
%
% Once the component parts have been created, they can be assembled
% to create the flame object.
%
fl = flame(gas, burner, f, s);
fl.setMaxJacAge(max_jacobian_age(1), max_jacobian_age(2));

% if the starting solution is to be read from a previously-saved
% solution, uncomment this line and edit the file name and solution id.
%restore(fl,'h2flame2.xml', 'energy')

fl.solve(loglevel, refine_grid);

%% Enable the energy equation
%
%  The energy equation will now be solved to compute the
%  temperature profile. We also tighten the grid refinement
%  criteria to get an accurate final solution.

f.enableEnergy;
fl.setRefineCriteria(2, 200.0, 0.05, 0.1);
fl.solve(1, 1);
fl.saveSoln('h2fl.xml', 'energy', ['solution with energy equation']);

%% Show statistics

fl.writeStats;
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);

%% Make plots

clf;
subplot(2, 2, 1);
plotSolution(fl, 'flow', 'T');
title('Temperature [K]');
subplot(2, 2, 2);
plotSolution(fl, 'flow', 'velocity');
title('Axial Velocity [m/s]');
subplot(2, 2, 3);
plotSolution(fl, 'flow', 'H2O');
title('H2O Mass Fraction');
subplot(2, 2, 4);
plotSolution(fl, 'flow', 'O2');
title('O2 Mass Fraction');

toc

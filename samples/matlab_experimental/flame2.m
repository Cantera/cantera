%% FLAME2 - An axisymmetric stagnation-point non-premixed flame
%
% This script simulates a stagnation-point ethane-air flame.
%
% Keywords: combustion, 1D flow, strained flame, diffusion flame, plotting

%% Initialization

clear all
close all

tic
help flame2

t0 = cputime; % record the starting time

%% Set parameter values

p = OneAtm; % pressure
tin = 300.0; % inlet temperature
mdot_o = 0.72; % air, kg/m^2/s
mdot_f = 0.24; % fuel, kg/m^2/s

rxnmech = 'gri30.yaml'; % reaction mechanism file
comp1 = 'O2:0.21, N2:0.78, AR:0.01'; % air composition
comp2 = 'C2H6:1'; % fuel composition

initial_grid = 0.02 * [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]; % m

tol_ss = {1.0e-5, 1.0e-13}; % {rtol atol} for steady-state
% problem
tol_ts = {1.0e-4, 1.0e-13}; % {rtol atol} for time stepping

loglevel = 1; % amount of diagnostic output (0 to 5)

refine_grid = 1; % 1 to enable refinement, 0 to disable

%% Create the gas object
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties

gas = Solution(rxnmech, 'gri30', 'mixture-averaged');

% set its state to that of the  fuel (arbitrary)
gas.TPX = {tin, p, comp2};

%% Create the flow object
%
f = AxisymmetricFlow(gas, 'flow');
f.P = p;
f.setupGrid(initial_grid);
f.setSteadyTolerances('default', tol_ss{:});
f.setTransientTolerances('default', tol_ts{:});

%% Create the air inlet
%
%  The temperature, mass flux, and composition (relative molar) may be
%  specified.

inlet_o = Inlet(gas, 'air_inlet');
inlet_o.T = tin;
inlet_o.massFlux = mdot_o;
inlet_o.setMoleFractions(comp1);

%% Create the fuel inlet
%
inlet_f = Inlet(gas, 'fuel_inlet');
inlet_f.T = tin;
inlet_f.massFlux = mdot_f;
inlet_f.setMoleFractions(comp2);

%% Create the flame object
%
% Once the component parts have been created, they can be assembled
% to create the flame object.

fl = flame(gas, inlet_o, f, inlet_f);

% if the starting solution is to be read from a previously-saved
% solution, uncomment this line and edit the file name and solution id.
%restore(fl,'h2flame2.xml', 'energy')

% solve with fixed temperature profile first
fl.solve(loglevel, refine_grid);

%% Enable the energy equation
%
%  The energy equation will now be solved to compute the
%  temperature profile. We also tighten the grid refinement
%  criteria to get an accurate final solution.

f.energyEnabled = true;
fl.setRefineCriteria(2, 200.0, 0.1, 0.1);
fl.solve(loglevel, refine_grid);

%% Show statistics

fl.writeStats;
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);

%% Make plots

figure(1);
subplot(2, 3, 1);
plotSolution(fl, 'flow', 'T');
title('Temperature [K]');
subplot(2, 3, 2);
plotSolution(fl, 'flow', 'C2H6');
title('C2H6 Mass Fraction');
subplot(2, 3, 3);
plotSolution(fl, 'flow', 'O2');
title('O2 Mass Fraction');
subplot(2, 3, 4);
plotSolution(fl, 'flow', 'CH');
title('CH Mass Fraction');
subplot(2, 3, 5);
plotSolution(fl, 'flow', 'spread_rate');
title('Radial Velocity / Radius [s^-1]');
subplot(2, 3, 6);
plotSolution(fl, 'flow', 'velocity');
title('Axial Velocity [m/s]');

toc

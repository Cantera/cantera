%% Axisymmetric stagnation-point non-premixed flame
%
% This script simulates a stagnation-point ethane-air flame.
%
% Requires: cantera >= 3.2.0
%
% .. tags:: Matlab, combustion, 1D flow, strained flame, diffusion flame, plotting

%% Initialization

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

width = 0.02;
nz = 11;

logLevel = 1; % amount of diagnostic output (0 to 5)
refineGrid = 1; % 1 to enable refinement, 0 to disable

%% Create the gas object
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties

gas = Solution(rxnmech, 'gri30', 'mixture-averaged');

% set its state to that of the fuel (arbitrary)
gas.TPX = {tin, p, comp2};

%% Create the flow object
%
f = AxisymmetricFlow(gas, 'flow');
f.P = p;
f.setupUniformGrid(nz, width, 0.0);

%% Create the air inlet
%
% The temperature, mass flux, and composition (relative molar) may be
% specified.

inlet_o = Inlet(gas, 'air_inlet');
inlet_o.T = tin;
inlet_o.massFlux = mdot_o;
inlet_o.X = comp1;

%% Create the fuel inlet
%
inlet_f = Inlet(gas, 'fuel_inlet');
inlet_f.T = tin;
inlet_f.massFlux = mdot_f;
inlet_f.X = comp2;

%% Create the flame object
%
% Once the component parts have been created, they can be assembled
% to create the flame object.

fl = flame(gas, inlet_o, f, inlet_f);

%%
% Set an initial temperature guess to ensure ignition

f.setProfile('T', [0.0, mdot_o/(mdot_o + mdot_f), 1.0], [tin, 2000, tin]);

% Solve with fixed temperature profile first
fl.solve(logLevel, refineGrid);

%% Enable the energy equation
%
% The energy equation will now be solved to compute the temperature profile. We also
% tighten the grid refinement criteria to get an accurate final solution.

f.energyEnabled = true;
f.setRefineCriteria(200.0, 0.1, 0.1);
fl.solve(logLevel, refineGrid);

%% Show statistics

fl.writeStats;
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);

%% Make plots

figure(1);
subplot(2, 3, 1);
plotSolution(f, 'T');
title('Temperature [K]');
subplot(2, 3, 2);
plotSolution(f, 'C2H6');
title('C2H6 Mass Fraction');
subplot(2, 3, 3);
plotSolution(f, 'O2');
title('O2 Mass Fraction');
subplot(2, 3, 4);
plotSolution(f, 'CH');
title('CH Mass Fraction');
subplot(2, 3, 5);
plotSolution(f, 'spreadRate');
title('Radial Velocity / Radius [s^-1]');
subplot(2, 3, 6);
plotSolution(f, 'velocity');
title('Axial Velocity [m/s]');

toc

%%
% Plotting Utility
% ----------------

function plotSolution(domain, component)
    % Utility for plotting a specific solution component
    z = domain.grid;
    x = domain.values(component);
    plot(z, x);
    xlabel('z (m)');
    ylabel(component);
end

%% CATCOMB - Catalytic combustion of a stagnation flow on a platinum surface
%
% This script solves a catalytic combustion problem. A stagnation flow
% is set up, with a gas inlet 10 cm from a platinum surface at 900
% K. The lean, premixed methane/air mixture enters at ~ 6 cm/s (0.06
% kg/m2/s), and burns catalytically on the platinum surface. Gas-phase
% chemistry is included too, and has some effect very near the
% surface.
%
% The catalytic combustion mechanism is from Deutschman et al., 26th
% Symp. (Intl.) on Combustion,1996 pp. 1747-1754
%
% Keywords: combustion, catalysis, 1D flow, surface chemistry

%% Initialization

help catcomb;

clear all
close all

t0 = cputime; % record the starting time

%% Set parameter values

p = oneatm; % pressure
tinlet = 300.0; % inlet temperature
tsurf = 900.0; % surface temperature
mdot = 0.06; % kg/m^2/s
transport = 'Mix'; % transport model

% Solve first for a hydrogen/air case for use as the initial estimate for
% the methane/air case.

% composition of the inlet premixed gas for the hydrogen/air case
comp1 = 'H2:0.05, O2:0.21, N2:0.78, AR:0.01';

% composition of the inlet premixed gas for the methane/air case
comp2 = 'CH4:0.095, O2:0.21, N2:0.78, AR:0.01';

% the initial grid, in meters. The inlet/surface separation is 10 cm.
initial_grid = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]; % m

% numerical parameters
tol_ss = {1.0e-8 1.0e-14}; % {rtol atol} for steady-state problem
tol_ts = {1.0e-4 1.0e-9}; % {rtol atol} for time stepping

loglevel = 1; % amount of diagnostic output
% (0 to 5)

refine_grid = 1; % 1 to enable refinement, 0 to
% disable

%% Create the gas object
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
% The gas phase will be taken from the definition of phase 'gas' in
% input file 'ptcombust.yaml', which is a stripped-down version of
% GRI-Mech 3.0.

gas = Solution('ptcombust.yaml', 'gas', transport);
gas.TPX = {tinlet, p, comp1};

%% Create the interface object
%
% This object will be used to evaluate all surface chemical production
% rates. It will be created from the interface definition 'Pt_surf'
% in input file 'ptcombust.yaml,' which implements the reaction
% mechanism of Deutschmann et al., 1995 for catalytic combustion on
% platinum.

surf_phase = Interface('ptcombust.yaml', 'Pt_surf', gas);
surf_phase.T = tsurf;

% integrate the coverage equations in time for 1 s, holding the gas
% composition fixed to generate a good starting estimate for the
% coverages.

surf_phase.advanceCoverages(1.0);

% The two objects we just created are independent of the problem
% type -- they are useful in zero-D simulations, 1-D simulations,
% etc. Now we turn to creating the objects that are specifically
% for 1-D simulations. These will be 'stacked' together to create
% the complete simulation.

%% Create the flow object
%
% The flow object is responsible for evaluating the 1D governing
% equations for the flow. We will initialize it with the gas
% object, and assign it the name 'flow'.

flow = AxisymmetricFlow(gas, 'flow');

% set some parameters for the flow
flow.P = p;
flow.setupGrid(initial_grid);
flow.setSteadyTolerances('default', tol_ss{:});
flow.setTransientTolerances('default', tol_ts{:});

%% create the inlet
%
%  The temperature, mass flux, and composition (relative molar) may be
%  specified. This object provides the inlet boundary conditions for
%  the flow equations.

inlt = Inlet('inlet');

% set the inlet parameters. Start with comp1 (hydrogen/air)
inlt.T = tinlet;
inlt.setMassFlowRate(mdot);
inlt.setMoleFractions(comp1);

%% create the surface
%
% This object provides the surface boundary conditions for the flow
% equations. By supplying object surface_phase as an argument, the
% coverage equations for its surface species will be added to the
% equation set, and used to compute the surface production rates of
% the gas-phase species.

surf = Surface('surface', surf_phase);
surf.T = tsurf;

%% create the stack
%
% Once the component parts have been created, they can be assembled
% to create the 1D simulation.

stack = Sim1D({inlt, flow, surf});

% set the initial profiles.
stack.setProfile(2, {'velocity', 'spread_rate', 'T'}, ...
                [0.0, 1.0 % z/zmax
                 0.06, 0.0 % u
                 0.0, 0.0 % V
                 tinlet, tsurf]); % T
names = gas.speciesNames;

for k = 1:gas.nSpecies
    y = inlt.massFraction(k);
    stack.setProfile(2, names{k}, [0, 1; y, y]);
end

stack

%setTimeStep(fl, 1.0e-5, [1, 3, 6, 12]);
%setMaxJacAge(fl, 4, 5);

%% Solution

% start with the energy equation on
flow.enableEnergy;

% disable the surface coverage equations, and turn off all gas and
% surface chemistry

surf.setCoverageEqs('off');
surf_phase.setMultiplier(0.0);
gas.setMultiplier(0.0);

% solve the problem, refining the grid if needed
stack.solve(1, refine_grid);

% now turn on the surface coverage equations, and turn the
% chemistry on slowly

surf.setCoverageEqs('on');

for iter = 1:6
    mult = 10.0^(iter - 6);
    surf_phase.setMultiplier(mult);
    gas.setMultiplier(mult);
    stack.solve(1, refine_grid);
end

% At this point, we should have the solution for the hydrogen/air
% problem. Now switch the inlet to the methane/air composition.
inlt.setMoleFractions(comp2);

% set more stringent grid refinement criteria
stack.setRefineCriteria(2, 100.0, 0.15, 0.2);

% solve the problem for the final time
stack.solve(loglevel, refine_grid);

% show the solution
stack

% save the solution
stack.saveSoln('catcomb.xml', 'energy', ['solution with energy equation']);

%% Show statistics

stack.writeStats;
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);

%% Make plots

clf;

subplot(3, 3, 1);
stack.plotSolution('flow', 'T');
title('Temperature [K]');

subplot(3, 3, 2);
stack.plotSolution('flow', 'velocity');
title('Axial Velocity [m/s]');

subplot(3, 3, 3);
stack.plotSolution('flow', 'spread_rate');
title('Radial Velocity / Radius [1/s]');

subplot(3, 3, 4);
stack.plotSolution('flow', 'CH4');
title('CH4 Mass Fraction');

subplot(3, 3, 5);
stack.plotSolution('flow', 'O2');
title('O2 Mass Fraction');

subplot(3, 3, 6);
stack.plotSolution('flow', 'CO');
title('CO Mass Fraction');

subplot(3, 3, 7);
stack.plotSolution('flow', 'CO2');
title('CO2 Mass Fraction');

subplot(3, 3, 8);
stack.plotSolution('flow', 'H2O');
title('H2O Mass Fraction');

subplot(3, 3, 9);
stack.plotSolution('flow', 'H2');
title('H2 Mass Fraction');

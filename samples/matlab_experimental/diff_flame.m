%% Opposed-flow diffusion flame
%
% This example uses the ``CounterFlowDiffusionFlame`` function to solve an
% opposed-flow diffusion flame for Ethane in Air. This example is the same
% as the :doc:`diffusion_flame.py <../python/onedim/diffusion_flame>`
% example without radiation.
%
% .. tags:: Matlab, combustion, 1D flow, diffusion flame, plotting

%% Initialization

tic % total running time of the script
help diff_flame

runtime = cputime; % Record the starting time

%% Parameter values of inlet streams

p = OneAtm; % Pressure
tin = 300.0; % Inlet temperature
mdot_o = 0.72; % Air mass flux, kg/m^2/s
mdot_f = 0.24; % Fuel mass flux, kg/m^2/s
transport = 'mixture-averaged'; % Transport model
% NOTE: Transport model needed if mechanism file does not have transport
% properties.

%% Set-up initial grid, loglevel, tolerances. Enable/Disable grid refinement

width = 0.2;
nz = 11;

tol_ss = {1.0e-5, 1.0e-9}; % {rtol atol} for steady-state problem
tol_ts = {1.0e-3, 1.0e-9}; % {rtol atol} for time stepping
logLevel = 1; % Amount of diagnostic output (0 to 5)
refineGrid = 1; % 1 to enable refinement, 0 to disable

%% Create the gas objects for the fuel and oxidizer streams
%
% These objects will be used to evaluate all thermodynamic, kinetic, and
% transport properties.

fuel = Solution('gri30.yaml', 'gri30', transport);
ox = Solution('gri30.yaml', 'gri30', transport);
oxcomp = 'O2:0.21, N2:0.78'; % Air composition
fuelcomp = 'C2H6:1'; % Fuel composition
% Set each gas mixture state with the corresponding composition.
fuel.TPX = {tin, p, fuelcomp};
ox.TPX = {tin, p, oxcomp};

%% Set-up the flow object
%
% For this problem, the ``AxisymmetricFlow`` model is needed. Set the state of
% the flow as the fuel gas object. This is arbitrary and is only used to
% initialize the flow object. Set the grid to the initial grid defined
% prior, same for the tolerances.

f = AxisymmetricFlow(fuel, 'flow');
f.P = p;
f.setupUniformGrid(nz, width, 0.0);
f.setSteadyTolerances(tol_ss{:}, 'default');
f.setTransientTolerances(tol_ts{:}, 'default');

%% Create the fuel and oxidizer inlet steams
%
% Specify the temperature, mass flux, and composition correspondingly.

% Set the oxidizer inlet.
inlet_o = Inlet(ox, 'air_inlet');
inlet_o.T = tin;
inlet_o.massFlux = mdot_o;
inlet_o.X = oxcomp;

% Set the fuel inlet.
inlet_f = Inlet(fuel, 'fuel_inlet');
inlet_f.T = tin;
inlet_f.massFlux = mdot_f;
inlet_f.X = fuelcomp;

%% Create the flame object
%
% Once the inlets have been created, they can be assembled
% to create the flame object. Function ``CounterFlorDiffusionFlame``
% (in ``Cantera/1D``) sets up the initial guess for the solution using a
% Burke-Schumann flame. The input parameters are: fuel inlet object, flow
% object, oxidizer inlet object, fuel gas object, oxidizer gas object, and
% the name of the oxidizer species as in character format.

fl = CounterFlowDiffusionFlame(inlet_f, f, inlet_o, fuel, ox, 'O2');

%% Solve with fixed temperature profile
%
% Grid refinement is turned off for this process in this example.
% To turn grid refinement on, change 0 to 1 for last input is solve function.

fl.solve(logLevel, 0);

%% Enable the energy equation
%
% The energy equation will now be solved to compute the temperature profile.
% We also tighten the grid refinement criteria to get an accurate final
% solution. The explanation of the ``setRefineCriteria`` function is located
% on cantera.org in the Matlab User's Guide and can be accessed by
% ``help setRefineCriteria``.

f.energyEnabled = true;
f.setRefineCriteria(200.0, 0.1, 0.2);
fl.solve(logLevel, refineGrid);

%% Show statistics of solution and elapsed time

fl.writeStats;
elapsed = cputime - runtime;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);

%% Plot results
% Make a single plot showing temperature and mass fraction of select
% species along axial distance from fuel inlet to air inlet.

figure(1);
subplot(2, 3, 1);
plotSolution(f, 'T');
title('Temperature [K]');
subplot(2, 3, 2);
plotSolution(f, 'velocity');
title('Axial Velocity [m/s]');
subplot(2, 3, 3);
plotSolution(f, 'O2');
title('O2 Mass Fraction');
subplot(2, 3, 4);
plotSolution(f, 'H2O');
title('H2O Mass Fraction');
subplot(2, 3, 5);
plotSolution(f, 'C2H6');
title('C2H6 Mass Fraction');
subplot(2, 3, 6);
plotSolution(f, 'CO2');
title('CO2 Mass Fraction');

toc

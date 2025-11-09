%% Opposed-flow diffusion flame
%
% This example uses the ``CounterFlowDiffusionFlame`` object to solve an
% opposed-flow diffusion flame for Ethane in Air. This example is the same
% as the :doc:`diffusion_flame.py <../python/onedim/diffusion_flame>`
% example without radiation.
%
% Requires: cantera >= 3.2.0
%
% .. tags:: Matlab, combustion, 1D flow, strained flame, diffusion flame, plotting

%%
% Problem Definition
% ------------------

tic  % total running time of the script
help diffusion_flame

runtime = cputime;  % Record the starting time

%%
% **Set parameter values**

p = ct.OneAtm;  % Pressure
tin = 300.0;  % Inlet temperature
mdot_o = 0.72;  % Air mass flux, kg/m^2/s
mdot_f = 0.24;  % Fuel mass flux, kg/m^2/s
transport = 'mixture-averaged';  % Transport model

%%
% Set-up initial grid, loglevel, tolerances. Enable/Disable grid refinement

width = 0.02;
nz = 11;

tol_ss = {1.0e-5, 1.0e-9};  % {rtol atol} for steady-state problem
tol_ts = {1.0e-3, 1.0e-9};  % {rtol atol} for time stepping
logLevel = 1;  % Amount of diagnostic output (0 to 5)
refineGrid = 1;  % 1 to enable refinement, 0 to disable

%%
% **Create the gas object for the fuel and oxidizer streams**
%
% This object will be used to evaluate all thermodynamic, kinetic, and
% transport properties.

gas = ct.Solution('gri30.yaml', 'gri30', transport);
oxcomp = 'O2:0.21, N2:0.78';  % Air composition
fuelcomp = 'C2H6:1';  % Fuel composition

%%
% **Set up the flow object**
%
% For this problem, the ``AxisymmetricFlow`` model is needed. Set the state of
% the flow as the fuel gas object. This is arbitrary and is only used to
% initialize the flow object. Set the grid to the initial grid defined
% prior, same for the tolerances.

flow = ct.oneD.AxisymmetricFlow(gas, 'flow');
flow.P = p;
flow.setupUniformGrid(nz, width, 0.0);
flow.setSteadyTolerances(tol_ss{:}, 'default');
flow.setTransientTolerances(tol_ts{:}, 'default');

%%
% **Create the fuel and oxidizer inlet steams**
%
% Specify the temperature, mass flux, and composition correspondingly.

% Set the fuel inlet.
gas.TPX = {tin, p, fuelcomp};
inlet_f = ct.oneD.Inlet(gas, 'fuel_inlet');
inlet_f.massFlux = mdot_f;

% Set the oxidizer inlet.
gas.TPX = {tin, p, oxcomp};
inlet_o = ct.oneD.Inlet(gas, 'air_inlet');
inlet_o.massFlux = mdot_o;

%%
% **Create the flame object**
%
% Once the inlets have been created, they can be assembled to create the flame object.
% The object ``CounterFlowDiffusionFlame`` sets up the initial guess for the solution
% using a Burke-Schumann flame. The input parameters are: fuel inlet object, flow
% object, oxidizer inlet object, and, optionally, the name of the oxidizer species.

fl = ct.oneD.CounterFlowDiffusionFlame(inlet_f, flow, inlet_o);

%%
% Solution
% --------
%
% **Solve with fixed temperature profile**
%
% Grid refinement is turned off for this process in this example.
% To turn grid refinement on, change 0 to 1 for last input is solve function.

fl.solve(logLevel, 0);

%%
% **Enable the energy equation**
%
% The energy equation will now be solved to compute the temperature profile.
% We also tighten the grid refinement criteria to get an accurate final
% solution. The explanation of the ``setRefineCriteria`` function is located
% on cantera.org in the Matlab User's Guide and can be accessed by
% ``help setRefineCriteria``.

flow.energyEnabled = true;
flow.setRefineCriteria(200.0, 0.1, 0.2);
fl.solve(logLevel, refineGrid);

%%
% **Show statistics and display results**

%%
% Show statistics of solution and elapsed time

fl.writeStats;
elapsed = cputime - runtime;
e = sprintf('Elapsed CPU time: %10.4g', elapsed);
disp(e);
toc

%%
% Make a single plot showing temperature and mass fraction of select
% species along axial distance from fuel inlet to air inlet.

figure(1);
subplot(2, 3, 1);
plotSolution(flow, 'T', 'Temperature [K]');
subplot(2, 3, 2);
plotSolution(flow, 'velocity', 'Axial Velocity [m/s]');
subplot(2, 3, 3);
plotSolution(flow, 'O2', 'O2 Mass Fraction');
subplot(2, 3, 4);
plotSolution(flow, 'H2O', 'H2O Mass Fraction');
subplot(2, 3, 5);
plotSolution(flow, 'C2H6', 'C2H6 Mass Fraction');
subplot(2, 3, 6);
plotSolution(flow, 'CO2', 'CO2 Mass Fraction');

%%
% Plotting Utility
% ----------------

function plotSolution(domain, component, titleStr)
    % Utility for plotting a specific solution component
    z = domain.grid;
    x = domain.values(component);
    plot(z, x);
    xlabel('z (m)');
    ylabel(component);
    title(titleStr);
end

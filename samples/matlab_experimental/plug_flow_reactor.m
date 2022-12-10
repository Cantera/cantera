% General_Plug_Flow_Reactor (PFR) - to solve PFR equations for reactors
%
%    This code snippet is to model a constant area and varying area
%    (converging and diverging) nozzle as Plug Flow Reactor with given
%    dimensions and an incoming gas. The pressure is not assumed to be
%    constant here, as opposed to the Python Version of the
%    Plug Flow Reactor model.
%
%    The reactor assumes that the flow follows the Ideal Gas Law.
%
%    The governing equations used in this code can be referenced at:
%    *S.R Turns, An Introduction to Combustion - Concepts and Applications,
%    McGraw Hill Education, India, 2012, 206-210.*
%
%    The current example is written for methane combustion, but can be readily
%    adapted for other chemistries.
%
%    Developed by Ashwin Kumar/Dr.Joseph Meadows (mgak@vt.edu/jwm84@vt.edu) on 3-June-2020
%    Research Assistant/Assistant Professor
%    Advanced Propulsion and Power Laboratory
%    Virginia Tech

%% Clear all variables, close all figures, clear the command line:

clear all
close all
cleanup
clc

tic
help PFR

%% Operation Parameters

% Temperature of gas, in K
T0 = 1473;
% Pressure of gas, in Pa
P0 = 4.47*101325;

% Equivalence Ratio
Phi = 0.2899;

% Import the gas phase, read out key species indices:
gas_calc = Solution('gri30.yaml', 'gri30');
ich4 = gas_calc.speciesIndex('CH4');
io2  = gas_calc.speciesIndex('O2');
in2  = gas_calc.speciesIndex('N2');
nsp = gas_calc.nSpecies;
x = zeros(nsp,1);

% Change the below values for different Phi values of methane Combustion
x(ich4,1) = Phi;
x(io2,1) = 2.0;
x(in2,1) = 7.52;

% Set the initial state and then equilibrate for a given enthalpy and pressure:
gas_calc.TPX = {T0, P0, x};
gas_calc.equilibrate('HP');

%% Calculation of properties along the reactor length

% The Dimensions and conditions of the reactor are given below

% Inlet Area, in m^2
A_in = 0.018;
% Exit Area, in m^2
A_out = 0.003;
% Length of the reactor, in m
L = 1.284*0.0254;
% The whole reactor is divided into n small reactors
n = 100;
% Mass flow rate into the reactor, in kg/s
mdot_calc = 1.125;

% Flag to indicate whether the area is converging, diverging, or constant:
% k = -1 makes the solver solve for converging area.
% k = +1 makes the solver solve for diverging area.
% k = 0 makes the solver solve for constant cross sectional area
if A_in > A_out
    k = -1;
elseif A_out > A_in
    k = 1;
else
    k = 0;
end

dAdx = abs(A_in-A_out)/L;
% The whole length of the reactor is divided into n small lengths
dx = L/n;

x_calc = 0:dx:L;
nsp = gas_calc.nSpecies;

% Initialize arrays for T, Y, and rho at each location:
T_calc = zeros(length(x_calc), 1);
Y_calc = zeros(length(x_calc), nsp);
rho_calc = zeros(length(x_calc), 1);

T_calc(1) = gas_calc.T;
Y_calc(1,:) = gas_calc.Y;
rho_calc(1) = gas_calc.D;

for i = 2:length(x_calc)

    %Solver location indicator
    fprintf('Solving reactor %d of %d\n', i, length(x_calc))

%--------------------------------------------------------------------------
%------The values of variables at previous location are given as initial---
%------values to the current iteration and the limits of the current-------
%--------------reactor and the gas entering it are being set---------------
%--------------------------------------------------------------------------
    inlet_soln(1) = rho_calc(i-1);
    inlet_soln(2) = T_calc(i-1);
    inlet_soln(3:nsp+2) = Y_calc(i-1, :);
    limits = [x_calc(i-1), x_calc(i)];
    gas_calc.TDY = {T_calc(i-1), rho_calc(i-1), Y_calc(i-1, :)};
    options = odeset('RelTol', 1.e-10, 'AbsTol', 1e-10,...
                     'InitialStep', 1e-8, 'NonNegative', 1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    % These values are passed onto the ode15s solver
    [~,y] = ode15s(@PFR_Solver, limits, inlet_soln, options, ...
                    gas_calc, mdot_calc, A_in, dAdx, k);

    T_calc(i) = y(end, 2);
    rho_calc(i) = y(end, 1);
    Y_calc(i,:) = y(end, 3:nsp+2);
end

A_calc = A_in + k.* x_calc * dAdx;
vx_calc = zeros(length(x_calc), 1);
R_calc = zeros(length(x_calc), 1);
M_calc = zeros(length(x_calc), 1);
P_calc = zeros(length(x_calc), 1);
for i=1:length(x_calc)
    % The gas is set to the solved property values at each location
    gas.TDY = {T_calc(i), rho_calc(i), Y_calc(i, :)};
    % Velocity is calculated from Mass flow rate, Area and Density
    vx_calc(i) = mdot_calc./ (A_calc(i) * rho_calc(i));
    % Specific Gas Constant
    R_calc(i) = gasconstant() / gas_calc.meanMolecularWeight;
    % Mach No. is calculated from local velocity and local speed of sound
    M_calc(i) = vx_calc(i) / gas_calc.soundspeed;
    % Pressure is calculated from density, temperature and gas constant
    P_calc(i) = rho_calc(i) * R_calc(i) * T_calc(i);
end

%% Plotting
plot(x_calc, M_calc)
xlabel('X-Position (m)')
ylabel('Mach No')
title('Mach No. Variation')
figure(2)
plot(x_calc, A_calc)
xlabel('X-Position (m)')
ylabel('Area (m^2)')
title('Reactor Profile')
figure(3)
plot(x_calc, T_calc)
xlabel('X-Position (m)')
ylabel('Temperature')
title('Temperature Variation')
figure(4)
plot(x_calc, rho_calc)
xlabel('X-Position (m)')
ylabel('Density (kg/m^3)')
title('Density Variation')
figure(5)
plot(x_calc, P_calc)
xlabel('X-Position (m)')
ylabel('Pressure (Pa)')
title('Pressure Variation')

toc

% General_Plug_Flow_Reactor(PFR) - to solve PFR equations for reactors
%
%    This code snippet is to model a constant area and varying area
%    (converging and diverging) nozzle as Plug Flow Reactor with given
%    dimensions and an incoming gas. The pressure is not assumed to be
%    constant here as opposed to the Python Version of the
%    Plug Flow Reactor modelling.
%
%    The governing equations used in this code can be referenced at:
%    *S.R Turns, An Introduction to Combustion - Concepts and Applications,
%    McGraw Hill Education, India, 2012, 206-210.*
%

%% Setting the Gas
clear all
close all
clc

%Temperature of gas
T0 = 1473;
%Pressure of gas
P0 = 4.47*101325;

%Equivalence Ratio
Phi = 0.2899;

gas_calc = Solution('gri30.yaml');
ich4 = speciesIndex(gas_calc,'CH4');
io2  = speciesIndex(gas_calc,'O2');
in2  = speciesIndex(gas_calc,'N2');
nsp = nSpecies(gas_calc);
x = zeros(nsp,1);
% Change the below values for different Phi values of methane Combustion
x(ich4,1) = Phi;
x(io2,1) = 2.0;
x(in2,1) = 7.52;
set(gas_calc,'T',T0,'P',P0,'MoleFractions',x);
gas_calc = equilibrate(gas_calc,'HP');

%% Calculation of properties along the reactor length
% The Dimensions and conditions of the reactor are given below

% Inlet Area
A_in = 0.018;
% Exit Area
A_out = 0.003;
% Length of the reactor
L = 1.284*0.0254;
% The whole reactor is divided into n small reactors
n = 100;
% Mass flow rate into the reactor in kg/s
mdot_calc = 1.125;
% k = -1 makes the solver solve for converging area.
% k = +1 makes the solver solve for diverging area.
% k = 0 makes the solver solve for constant cross sectional area
if A_in>A_out
    k = -1;
elseif A_out>A_in
    k = 1;
else
    k = 0;
end

dAdx = abs(A_in-A_out)/L;
% The whole length of the reactor is divided into n small lengths
dx = L/n;

x_calc = 0:dx:L;
nsp = nSpecies(gas_calc);

T_calc = zeros(length(x_calc),1);
Y_calc = zeros(length(x_calc), nsp);
rho_calc = zeros(length(x_calc),1);

T_calc(1) = temperature(gas_calc);
Y_calc(1,:) = massFractions(gas_calc);
rho_calc(1) = density(gas_calc);

for i = 2:length(x_calc)

    %Solver location indicator
    fprintf('Solving reactor %d of %d\n', i, length(x_calc))

%--------------------------------------------------------------------------
%------The values of variables at previous location are given as initial---
%------values to the current iteration and the limits of the current-------
%--------------reactor and the gas entering it are being set---------------
%--------------------------------------------------------------------------
    init(1) = rho_calc(i-1);
    init(2) = T_calc(i-1);
    init(3:nsp+2) = Y_calc(i-1,:);
    limits = [x_calc(i-1),x_calc(i)];
    set(gas_calc,'T',T_calc(i-1),'Density',rho_calc(i-1),'MoleFractions',Y_calc(i-1,:));
    options = odeset('RelTol',1.e-10,'AbsTol',1e-10,'InitialStep',1e-8,'NonNegative',1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    % These values are passed onto the ode15s solver
    [~,y] = ode15s(@PFR_solver,limits,init,options,gas_calc,mdot_calc,A_in,dAdx,k);

    T_calc(i) = y(end,2);
    rho_calc(i) = y(end,1);
    Y_calc(i,:) = y(end,3:nsp+2);
end

A_calc = A_in+k.*x_calc*dAdx;
vx_calc = zeros(length(x_calc),1);
R_calc = zeros(length(x_calc),1);
M_calc = zeros(length(x_calc),1);
P_calc = zeros(length(x_calc),1);
for i=1:length(x_calc)
    % The gas is set to the solved property values at each location
    set(gas_calc,'Temperature',T_calc(i),'Density',rho_calc(i),'MassFractions',Y_calc(i,:));
    % Velocity is calculated from Mass flow rate, Area and Density
    vx_calc(i) = mdot_calc./(A_calc(i)*rho_calc(i));
    % Specific Gas Constant
    R_calc(i) = gasconstant()/meanMolecularWeight(gas_calc);
    % Mach No. is calculated from local velocity and local speed of sound
    M_calc(i) = vx_calc(i)/soundspeed(gas_calc);
    % Pressure is calculated from density, temeprature and gas constant
    P_calc(i) = rho_calc(i)*R_calc(i)*T_calc(i);
end

%% Plotting
plot(x_calc,M_calc)
xlabel('X-Position (m)')
ylabel('Mach No')
title('Mach No. Variation')
figure(2)
plot(x_calc,A_calc)
xlabel('X-Position (m)')
ylabel('Area (m^2)')
title('Reactor Profile')
figure(3)
plot(x_calc,T_calc)
xlabel('X-Position (m)')
ylabel('Temperature')
title('Temperature Variation')
figure(4)
plot(x_calc,rho_calc)
xlabel('X-Position (m)')
ylabel('Density (kg/m^3)')
title('Density Variation')
figure(5)
plot(x_calc,P_calc)
xlabel('X-Position (m)')
ylabel('Pressure (Pa)')
title('Pressure Variation')

%% Freely-propagating flame
%
% MATLAB demo program to compute flame speeds using GRI-Mech.
%
% .. tags:: Matlab, combustion, 1D flow, premixed flame, flame speed, saving output

%% Initialization

tic
help flamespeed

t0 = cputime;  % record the starting time

%% Set parameter values

phi = 1.;  % equivalence ratio
Tin = 300.0;  % burner temperature (K)
pressure = OneAtm;  % pressure

uIn = 0.3;  % initial guess (m/s)
lz = 0.1;  % length of simulated flow domain
nz = 6;  % initial number of grid points

logLevel = 1;
refineGrid = true;  % turn on grid refinement

%% Create the gas object
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties

gas = Solution('gri30.yaml', 'gri30', 'mixture-averaged');
nsp = gas.nSpecies;

gas.setEquivalenceRatio(phi, 'CH4', 'O2:0.21,N2:0.79');
gas.TP = {Tin, pressure};
rhoIn = gas.massDensity;
Xin = gas.X;
Yin = gas.Y;

gas.equilibrate('HP');
rhoOut = gas.massDensity;
Yout = gas.Y;
Tad = gas.T;
disp(sprintf("phi = %1.1f, Tad = %1.1f\n", phi, Tad));

%% Create domains required for a flame simulation
%

% Create the flow domain
flame = FreeFlow(gas);
flame.setupUniformGrid(nz, lz, 0.);
flame.setSteadyTolerances(1e-5, 1e-11);

% Create the inlet
inlt = Inlet(gas);
mdot = uIn * rhoIn;

% Create the outlet
outlt = Outlet(gas);
uOut = uIn * rhoIn / rhoOut;

%% Create the stack and insert the domains
%

stack = Sim1D({inlt, flame, outlt});

% Supply initial guess

locs = [0.0, 0.3, 0.7, 1.0];
flame.setProfile('velocity', locs, [uIn, uIn, uOut, uOut]);
flame.setProfile('T', locs, [Tin, Tin, Tad, Tad]);

names = gas.speciesNames;
for i = 1:gas.nSpecies
    flame.setProfile(names{i}, locs, [Yin(i), Yin(i), Yout(i), Yout(i)]);
end
flame.P = pressure;

inlt.massFlux = mdot;
inlt.T = Tin;
inlt.X = Xin;

stack.show();

ratio = 10.0;
slope = 0.08;
curve = 0.1;

flame.setRefineCriteria(ratio,slope,curve);

fprintf("Profile used for initial guess:\n\n")
flame.info

%% Save the initial guess to a container file
%

if ctUsesHDF5()
    % Cantera is compiled with native HDF5 support
    fileName = 'flamespeed.h5';
else
    fileName = 'flamespeed.yaml';
end
stack.save(fileName, 'initial-guess', 'Initial guess', true);

%% Solve freely propagating flame
%

% Linearly interpolate to find location where this temperature would
% exist. The temperature at this location will then be fixed for
% remainder of calculation.
stack.setFixedTemperature(0.5 * (Tin + Tad));

flame.energyEnabled = true;
stack.solve(logLevel, refineGrid);
uVec = flame.values('velocity');
disp(sprintf("Flame speed with mixture-averaged transport: %1.4f m/s\n", uVec(1)));
stack.save(fileName, "mix", "Solution with mixture-averaged transport", true);

% now switch to multicomponent transport
flame.transportModel = 'multicomponent';
stack.solve(logLevel, refineGrid);
uVec = flame.values('velocity');
disp(sprintf("Flame speed with multicomponent transport: %1.4f m/s\n", uVec(1)));
stack.save(fileName, "multi", "Solution with multicomponent transport", true);

% now enable Soret diffusion
flame.soretEnabled = true;
stack.solve(logLevel, refineGrid);
uVec = flame.values('velocity');
disp(sprintf("Flame speed with multicomponent transport + Soret: %1.4f m/s\n", uVec(1)));
stack.save(fileName, "soret", "Solution with multicomponent transport and Soret", true);

fprintf("Profile after final solution step:\n\n")
flame.info

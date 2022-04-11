% Catalytic combustion of a stagnation flow on a platinum surface
%
% This script solves a catalytic combustion problem. A stagnation flow
% is set up, with a gas inlet 10 cm from a platinum surface at 900
% K. The lean, premixed methane/air mixture enters at ~ 6 cm/s (0.06
% kg/m2/s), and burns catalytically on the platinum surface. Gas-phase
% chemistry is included too, and has some effect very near the
% surface.
%
% The catalytic combustion mechanism is from Deutschmann et al., 26th
% Symp. (Intl.) on Combustion,1996 pp. 1747-1754
%
% Keywords: combustion, catalysis, 1D flow, surface chemistry

help catcomb;

clear all;
cleanup;

t0 = cputime;  % record the starting time

% Parameter values are collected here to make it easier to modify
% them

p          =   oneatm;              % pressure
tinlet     =   300.0;               % inlet temperature
tsurf      =   900.0;               % surface temperature
mdot       =   0.06;                % kg/m^2/s
transport  =  'Mix';                % transport model

% We will solve first for a hydrogen/air case to
% use as the initial estimate for the methane/air case

% composition of the inlet premixed gas for the hydrogen/air case
comp1       =  'H2:0.05, O2:0.21, N2:0.78, AR:0.01';

% composition of the inlet premixed gas for the methane/air case
comp2       =  'CH4:0.095, O2:0.21, N2:0.78, AR:0.01';

% the initial grid, in meters. The inlet/surface separation is 10 cm.
initial_grid = [0.0 0.02 0.04 0.06 0.08 0.1];  % m

% numerical parameters
tol_ss    = [1.0e-8 1.0e-14];       % [rtol atol] for steady-state problem
tol_ts    = [1.0e-4 1.0e-9];       % [rtol atol] for time stepping

loglevel  = 1;                      % amount of diagnostic output
                                    % (0 to 5)

refine_grid = 1;                    % 1 to enable refinement, 0 to
                                    % disable

%%%%%%%%%%%%%%%  end of parameter list  %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% create the gas object %%%%%%%%%%%%%%%%%%%%%%%%
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
% The gas phase will be taken from the definition of phase 'gas' in
% input file 'ptcombust.yaml,' which is a stripped-down version of
% GRI-Mech 3.0.
gas = Solution('ptcombust.yaml','gas',transport);
set(gas,'T',tinlet,'P',p,'X',comp1);

%%%%%%%%%%%%%%%% create the interface object %%%%%%%%%%%%%%%%%%
%
% This object will be used to evaluate all surface chemical production
% rates. It will be created from the interface definition 'Pt_surf'
% in input file 'ptcombust.yaml,' which implements the reaction
% mechanism of Deutschmann et al., 1995 for catalytic combustion on
% platinum.
%
surf_phase = importInterface('ptcombust.yaml','Pt_surf',gas);
setTemperature(surf_phase, tsurf);

% integrate the coverage equations in time for 1 s, holding the gas
% composition fixed to generate a good starting estimate for the
% coverages.
advanceCoverages(surf_phase, 1.0);

% The two objects we just created are independent of the problem
% type -- they are useful in zero-D simulations, 1-D simulations,
% etc. Now we turn to creating the objects that are specifically
% for 1-D simulations. These will be 'stacked' together to create
% the complete simulation.

%%%%%%%%%%%%%%%% create the flow object %%%%%%%%%%%%%%%%%%%%%%%
%
% The flow object is responsible for evaluating the 1D governing
% equations for the flow. We will initialize it with the gas
% object, and assign it the name 'flow'.
%
flow = AxisymmetricFlow(gas, 'flow');

% set some parameters for the flow
set(flow, 'P', p, 'grid', initial_grid, 'tol', tol_ss, 'tol-time', tol_ts);

%%%%%%%%%%%%%%% create the inlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The temperature, mass flux, and composition (relative molar) may be
%  specified. This object provides the inlet boundary conditions for
%  the flow equations.
%
inlt = Inlet('inlet');

% set the inlet parameters. Start with comp1 (hydrogen/air)
set(inlt, 'T', tinlet, 'MassFlux', mdot, 'X', comp1);

%%%%%%%%%%%%%% create the surface %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This object provides the surface boundary conditions for the flow
% equations. By supplying object surface_phase as an argument, the
% coverage equations for its surface species will be added to the
% equation set, and used to compute the surface production rates of
% the gas-phase species.
%
surf = Surface('surface', surf_phase);
setTemperature(surf,tsurf);

%%%%%%%%%%%%% create the stack  %%%%%%%%%%%%
%
% Once the component parts have been created, they can be assembled
% to create the 1D simulation.
%
sim1D = Stack([inlt, flow, surf]);

% set the initial profiles.
setProfile(sim1D, 2, {'velocity', 'spread_rate', 'T'}, ...
           [0.0            1.0       % z/zmax
            0.06           0.0       % velocity (u)
            0.0            0.0       % spread rate (V)
            tinlet         tsurf]);  % T
names = speciesNames(gas);
for k = 1:nSpecies(gas)
  y = massFraction(inlt, k);
  setProfile(sim1D, 2, names{k}, [0    1; y y]);
end
sim1D

%setTimeStep(fl, 1.0e-5, [1, 3, 6, 12]);
%setMaxJacAge(fl, 4, 5);

%%%%%%%%%%%%% solution %%%%%%%%%%%%%%%%%%%%

% start with the energy equation on
enableEnergy(flow);

% disable the surface coverage equations, and turn off all gas and
% surface chemistry
setCoverageEqs(surf, 'off');
setMultiplier(surf_phase, 0.0);
setMultiplier(gas, 0.0);

% solve the problem, refining the grid if needed
solve(sim1D, 1, refine_grid);

% now turn on the surface coverage equations, and turn the
% chemistry on slowly
setCoverageEqs(surf, 'on');
for iter=1:6
  mult = 10.0^(iter - 6);
  setMultiplier(surf_phase, mult);
  setMultiplier(gas, mult);
  solve(sim1D, 1, refine_grid);
end

% At this point, we should have the solution for the hydrogen/air
% problem. Now switch the inlet to the methane/air composition.
set(inlt,'X',comp2);

% set more stringent grid refinement criteria
setRefineCriteria(sim1D, 2, 100.0, 0.15, 0.2);

% solve the problem for the final time
solve(sim1D, loglevel, refine_grid);

% show the solution
sim1D

% save the solution
saveSoln(sim1D,'catcomb.yaml','energy',['solution with energy equation']);

%%%%%%%%%% show statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeStats(sim1D);
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g',elapsed);
disp(e);

%%%%%%%%%% make plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
subplot(3,3,1);
plotSolution(sim1D, 'flow', 'T');
title('Temperature [K]');
subplot(3,3,2);
plotSolution(sim1D, 'flow', 'u');
title('Axial Velocity [m/s]');
subplot(3,3,3);
plotSolution(sim1D, 'flow', 'V');
title('Radial Velocity / Radius [1/s]');
subplot(3,3,4);
plotSolution(sim1D, 'flow', 'CH4');
title('CH4 Mass Fraction');
subplot(3,3,5);
plotSolution(sim1D, 'flow', 'O2');
title('O2 Mass Fraction');
subplot(3,3,6);
plotSolution(sim1D, 'flow', 'CO');
title('CO Mass Fraction');
subplot(3,3,7);
plotSolution(sim1D, 'flow', 'CO2');
title('CO2 Mass Fraction');
subplot(3,3,8);
plotSolution(sim1D, 'flow', 'H2O');
title('H2O Mass Fraction');
subplot(3,3,9);
plotSolution(sim1D, 'flow', 'H2');
title('H2 Mass Fraction');

% FLAME1 - A burner-stabilized flat flame
%
%    This script simulates a burner-stablized lean hydrogen-oxygen flame
%    at low pressure.
%
% Keywords: combustion, 1D flow, burner-stabilized flame, plotting

help flame1

t0 = cputime;  % record the starting time

% parameter values
p          =   0.05*oneatm;         % pressure
tburner    =   373.0;               % burner temperature
mdot       =   0.06;                % kg/m^2/s

rxnmech    =  'h2o2.yaml';           % reaction mechanism file
comp       =  'H2:1.8, O2:1, AR:7'; % premixed gas composition

initial_grid = [0.0 0.02 0.04 0.06 0.08 0.1 ...
                0.15 0.2 0.4 0.49 0.5];  % m

tol_ss    = [1.0e-5 1.0e-13];       % [rtol atol] for steady-state
                                    % problem
tol_ts    = [1.0e-4 1.0e-9];        % [rtol atol] for time stepping

loglevel  = 1;                      % amount of diagnostic output (0
                                    % to 5)

refine_grid = 1;                    % 1 to enable refinement, 0 to
                                    % disable
max_jacobian_age = [5, 10];

%%%%%%%%%%%%%%%% create the gas object %%%%%%%%%%%%%%%%%%%%%%%%
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
gas = Solution(rxnmech, 'ohmech', 'mixture-averaged');

% set its state to that of the unburned gas at the burner
set(gas,'T', tburner, 'P', p, 'X', comp);

%%%%%%%%%%%%%%%% create the flow object %%%%%%%%%%%%%%%%%%%%%%%

f = AxisymmetricFlow(gas,'flow');
set(f, 'P', p, 'grid', initial_grid);
set(f, 'tol', tol_ss, 'tol-time', tol_ts);

%%%%%%%%%%%%%%% create the burner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The burner is an Inlet object. The temperature, mass flux,
%  and composition (relative molar) may be specified.
%
burner = Inlet('burner');
set(burner, 'T', tburner, 'MassFlux', mdot, 'X', comp);

%%%%%%%%%%%%%% create the outlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The type of flame is determined by the object that terminates
%  the domain. An Outlet object imposes zero gradient boundary
%  conditions for the temperature and mass fractions, and zero
%  radial velocity and radial pressure gradient.
%
s = Outlet('out');

%%%%%%%%%%%%% create the flame object  %%%%%%%%%%%%
%
% Once the component parts have been created, they can be assembled
% to create the flame object.
%
fl = flame(gas, burner, f, s);
setMaxJacAge(fl, max_jacobian_age(1),  max_jacobian_age(2));

% if the starting solution is to be read from a previously-saved
% solution, uncomment this line and edit the file name and solution id.
%restore(fl,'h2flame.yaml', 'energy')

solve(fl, loglevel, refine_grid);

%%%%%%%%%%%% enable the energy equation %%%%%%%%%%%%%%%%%%%%%
%
%  The energy equation will now be solved to compute the
%  temperature profile. We also tighten the grid refinement
%  criteria to get an accurate final solution.
%

enableEnergy(f);
setRefineCriteria(fl, 2, 200.0, 0.05, 0.1);
solve(fl, 1, 1);
saveSoln(fl,'h2flame.yaml','energy',['solution with energy equation']);

%%%%%%%%%% show statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeStats(fl);
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g',elapsed);
disp(e);

%%%%%%%%%% make plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
subplot(2,2,1);
plotSolution(fl, 'flow', 'T');
title('Temperature [K]');
subplot(2,2,2);
plotSolution(fl, 'flow', 'u');
title('Axial Velocity [m/s]');
subplot(2,2,3);
plotSolution(fl, 'flow', 'H2O');
title('H2O Mass Fraction');
subplot(2,2,4);
plotSolution(fl, 'flow', 'O2');
title('O2 Mass Fraction');

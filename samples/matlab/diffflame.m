% DIFFFLAME - An opposed-flow diffusion flame.
% This example uses the CounterFlowDiffusionFlame function to solve an
% opposed-flow diffusion flame for Ethane in Air. This example is the same
% as the diffusion_flame.py example without radiation.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

runtime = cputime;  % Record the starting time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter values of inlet streams
%

p          =   oneatm;              % Pressure
tin        =   300.0;               % Inlet temperature
mdot_o     =   0.72;                % Air mass flux, kg/m^2/s
mdot_f     =   0.24;                % Fuel mass flux, kg/m^2/s
transport  =  'Mix';                % Transport model
% NOTE: Transport model needed if mechanism file does not have transport
% properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set-up initial grid, loglevel, tolerances. Enable/Disable grid
% refinement.
%

initial_grid = 0.02*[0.0 0.2 0.4 0.6 0.8 1.0];  % Units: m
tol_ss    = [1.0e-5 1.0e-9];        % [rtol atol] for steady-state problem
tol_ts    = [1.0e-3 1.0e-9];        % [rtol atol] for time stepping
loglevel  = 1;                      % Amount of diagnostic output (0 to 5)
refine_grid = 1;                    % 1 to enable refinement, 0 to disable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the gas objects for the fuel and oxidizer streams. These objects
% will be used to evaluate all thermodynamic, kinetic, and transport
% properties.
%

fuel = GRI30(transport);
ox = GRI30(transport);
oxcomp     =  'O2:0.21, N2:0.78';   % Air composition
fuelcomp   =  'C2H6:1';             % Fuel composition
% Set each gas mixture state with the corresponding composition.
set(fuel,'T', tin, 'P', p, 'X', fuelcomp);
set(ox,'T',tin,'P',p,'X', oxcomp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set-up the flow object. For this problem, the AxisymmetricFlow model is
% needed. Set the state of the flow as the fuel gas object. This is
% arbitrary and is only used to initialize the flow object. Set the grid to
% the initial grid defined prior, same for the tolerances.
%

f = AxisymmetricFlow(fuel,'flow');
set(f, 'P', p, 'grid', initial_grid);
set(f, 'tol', tol_ss, 'tol-time', tol_ts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the fuel and oxidizer inlet steams. Specify the temperature, mass
% flux, and composition correspondingly.
%

% Set the oxidizer inlet.
inlet_o = Inlet('air_inlet');
set(inlet_o, 'T', tin, 'MassFlux', mdot_o, 'X', oxcomp);
%
% Set the fuel inlet.
inlet_f = Inlet('fuel_inlet');
set(inlet_f, 'T', tin, 'MassFlux', mdot_f, 'X', fuelcomp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Once the inlets have been created, they can be assembled
% to create the flame object. Function CounterFlorDiffusionFlame
% (in Cantera/1D) sets up the initial guess for the solution using a
% Burke-Schumann flame. The input parameters are: fuel inlet object, flow
% object, oxidizer inlet object, fuel gas object, oxidizer gas object, and
% the name of the oxidizer species as in character format.
%

fl = CounterFlowDiffusionFlame(inlet_f, f, inlet_o, fuel, ox, 'O2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve with fixed temperature profile first. Grid refinement is turned off
% for this process in this example. To turn grid refinement on, change 0 to
% 1 for last input is solve function.
%

solve(fl, loglevel, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enable the energy equation. The energy equation will now be solved to
% compute the temperature profile. We also tighten the grid refinement
% criteria to get an accurate final solution. The explanation of the
% setRefineCriteria function is located on cantera.org in the Matlab User's
% Guide and can be accessed by help setRefineCriteria
%

enableEnergy(f);
setRefineCriteria(fl, 2, 200.0, 0.1, 0.2);
solve(fl, loglevel, refine_grid);
saveSoln(fl,'c2h6.xml','energy',['solution with energy equation']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show statistics of solution and elapsed time.
%

writeStats(fl);
elapsed = cputime - runtime;
e = sprintf('Elapsed CPU time: %10.4g',elapsed);
disp(e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a single plot showing temperature and mass fraction of select
% species along axial distance from fuel inlet to air inlet.
%

z = grid(fl, 'flow');                    % Get grid points of flow
spec = speciesNames(fuel);               % Get species names in gas
T = solution(fl, 'flow', 'T');           % Get temperature solution
for i = 1:length(spec)
    % Get mass fraction of all species from solution
    y(i,:) = solution(fl, 'flow', spec{i});
end
j = speciesIndex(fuel, 'O2');            % Get index of O2 in gas object
k = speciesIndex(fuel, 'H2O');           % Get index of H2O in gas object
l = speciesIndex(fuel, 'C2H6');          % Get index of C2H6 in gas object
m = speciesIndex(fuel, 'CO2');           % Get index of CO2 in gas object
clf;
yyaxis left
plot(z,T)
xlabel('z (m)');
ylabel('Temperature (K)');
yyaxis right
plot(z,y(j,:),'r',z,y(k,:),'g',z,y(l,:),'m',z,y(m,:),'b');
ylabel('Mass Fraction');
legend('T','O2','H2O','C2H6','CO2');

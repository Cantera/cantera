%% Zero-dimensional reactor with surface chemistry
%
% This example illustrates how to use class ``Reactor`` for zero-dimensional
% simulations including both homogeneous and heterogeneous chemistry.
%
% Requires: cantera >= 3.2.0
%
% .. tags:: Matlab, catalysis, combustion, reactor network, plotting

%%
% Problem Definition
% ------------------

tic
help surf_reactor

% The surface reaction mechanism describes catalytic combustion of
% methane on platinum, and is from Deutschmann et al., 26th
% Symp. (Intl.) on Combustion,1996, pp. 1747-1754
surf = ct.Interface('ptcombust.yaml', 'Pt_surf');
gas = surf.adjacent('gas');

%%
% Set the initial conditions

t = 870.0;
gas.TPX = {t, ct.OneAtm, 'CH4:0.01, O2:0.21, N2:0.78'};

surf.TP = {t, surf.P};

nsp = gas.nSpecies;
nSurfSp = surf.nSpecies;

%%
% Create a reactor, and insert the gas
r = ct.IdealGasReactor(gas, 'reactor');
r.V = 1.0e-6;

%%
% Create a reservoir to represent the environment
a = ct.Solution('air.yaml', 'air', 'none');
a.TP = {t, ct.OneAtm};
env = ct.Reservoir(a);

%%
% Define a wall between the reactor and the environment and
% make it flexible, so that the pressure in the reactor is held
% at the environment pressure.
w = ct.Wall(r, env);

A = 1e-4; % Wall area

%%
% Add a reacting surface, with an area matching that of the wall
rsurf = ct.ReactorSurface(surf, r, 'surface');
rsurf.area = A;
rphase = rsurf.phase  % output needs to use phase owned by reactor

%%
% Set the wall area and heat transfer coefficient.
w.area = A;
w.heatTransferCoeff = 1.0e1; % W/m2/K

%%
% Set expansion rate parameter. dV/dt = KA(P_1 - P_2)
w.expansionRateCoeff = 1.0;

network = ct.ReactorNet({r});

%%
% Solution
% --------

nSteps = 100;
p0 = r.P;
names = {'CH4', 'CO', 'CO2', 'H2O'};
x = zeros([nSteps, 4]);
tim = zeros(nSteps, 1);
temp = zeros(nSteps, 1);
pres = zeros(nSteps, 1);
cov = zeros([nSteps nSurfSp]);
t = 0;
dt = 0.1;
t0 = cputime;

for n = 1:nSteps
    t = t + dt;
    network.advance(t);
    tim(n) = t;
    temp(n) = r.T;
    pres(n) = r.P - p0;
    cov(n, :) = rphase.X';
    x(n, :) = r.phase.moleFraction(names);
end

disp(['CPU time = ' num2str(cputime - t0)]);

%%
% Plot results

clf;
subplot(2, 2, 1);
plot(tim, temp);
xlabel('Time (s)');
ylabel('Temperature (K)');
subplot(2, 2, 2);
plot(tim, pres);
axis([0 5 -0.1 0.1]);
xlabel('Time (s)');
ylabel('Delta Pressure (Pa)');
subplot(2, 2, 3);
semilogy(tim, cov);
xlabel('Time (s)');
ylabel('Coverages');
legend(rphase.speciesNames);
subplot(2, 2, 4);
plot(tim, x);
xlabel('Time (s)');
ylabel('Mole Fractions');
legend(names);

toc

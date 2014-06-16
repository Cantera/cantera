% SURFREACTOR Zero-dimensional reactor with surface chemistry
% 
%    This example illustrates how to use class 'Reactor' for
%    zero-dimensional simulations including both homogeneous and
%    heterogeneous chemistry. 

help surfreactor

t = 870.0;
gas = importPhase('ptcombust.cti','gas');

% set the initial conditions
set(gas,'T',t,'P',oneatm,'X','CH4:0.01, O2:0.21, N2:0.78');

% The surface reaction mechanism describes catalytic combustion of
% methane on platinum, and is from Deutschman et al., 26th
% Symp. (Intl.) on Combustion,1996, pp. 1747-1754
surf = importInterface('ptcombust.cti','Pt_surf', gas);
setTemperature(surf, t);

nsp = nSpecies(gas);

% create a reactor, and insert the gas
r = IdealGasReactor(gas);
setInitialVolume(r, 1.0e-6)

% create a reservoir to represent the environment
a = IdealGasMix('air.cti');
set(a,'T',t,'P',oneatm);
env = Reservoir(a);

% Define a wall between the reactor and the environment and
% make it flexible, so that the pressure in the reactor is held
% at the environment pressure.
w = Wall;
install(w,r,env);

% set the surface mechanism on the left side of the wall (facing
% reactor 'r' to 'surf'. No surface mechanism will be installed on
% the air side.
setKinetics(w, surf, 0);

% set the wall area and heat transfer coefficient.
setArea(w, 1.0e-4);
setHeatTransferCoeff(w,1.0e1);  % W/m2/K

% set expansion rate parameter. dV/dt = KA(P_1 - P_2)
setExpansionRateCoeff(w, 1.0);

network = ReactorNet({r});
% setTolerances(network, 1.0e-8, 1.0e-12);

t = 0;
dt = 0.1;
t0 = cputime;
p0 = pressure(r);
names = {'CH4','CO','CO2','H2O'};
x = zeros([100 4]);
for n = 1:100
  t = t + dt;
  advance(network, t);
  tim(n) = t;
  temp(n) = temperature(r);  
  pres(n) = pressure(r) - p0;
  cov(n,:) = coverages(surf)'; 
  x(n,:) = moleFraction(gas,names);
end
disp(['CPU time = ' num2str(cputime - t0)]);

clf;
subplot(2,2,1);
plot(tim,temp);
xlabel('Time (s)');
ylabel('Temperature (K)');
subplot(2,2,2);
plot(tim,pres);
axis([0 5 -0.1 0.1]);
xlabel('Time (s)');
ylabel('Delta Pressure (Pa)');
subplot(2,2,3);
semilogy(tim,cov);
xlabel('Time (s)');
ylabel('Coverages');
legend(speciesNames(surf));
subplot(2,2,4);
plot(tim,x);
xlabel('Time (s)');
ylabel('Mole Fractions');
legend(names);
clear all
cleanup

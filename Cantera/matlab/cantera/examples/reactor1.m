function reactor1(g)
% REACTOR1 Zero-dimensional kinetics: adiabatic, constant pressure.
% 
%    This example illustrates how to use class 'Reactor' for
%    zero-dimensional kinetics simulations. Here the parameters are
%    set so that the reactor is adiabatic and very close to constant
%    pressure.
%

help reactor1

if nargin == 1 & isa(g,'solution')
   gas = g;
else
   gas = GRI30;
end

nsp = nSpecies(gas);

% set the initial conditions
set(gas,'T',1001.0,'P',oneatm,'X','H2:2,O2:1,N2:4');

% create a reactor, and insert the gas
r = Reactor;
insert(r, gas);

% create a reservoir to represent the environment
env = Reservoir;
a = IdealGasMix('air.xml');
insert(env, a);

% Define a wall between the reactor and the environment, and
% make it flexible, so that the pressure in the reactor is held
% at the environment pressure.
w = Wall;
install(w,r,env);

% set expansion parameter. dV/dt = K(P_1 - P_2)
setExpansionRateCoeff(w, 1.0e6);

% set wall area
setArea(w, 1.0);

t = 0;
dt = 1.0e-5;
t0 = cputime;
for n = 1:100
  t = t + dt;
  advance(r, t);
  disp([time(r) temperature(r)]);  
end
disp(['CPU time = ' num2str(cputime - t0)]);

clear all

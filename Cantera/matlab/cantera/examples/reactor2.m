function reactor2(g)
% REACTOR2 Zero-dimensional kinetics: adiabatic, constant volume.
% 
%    This example illustrates how to use class 'Reactor' for
%    zero-dimensional kinetics simulations. Here the parameters are
%    set so that the reactor is adiabatic and constant volume.
%

help reactor2

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

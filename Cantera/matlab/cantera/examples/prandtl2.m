function prandtl2(g)
% PRANDTL2  Prandlt number for an equilibrium H/O gas mixture.
%
%    This example does the same thing as prandtl1, but using
%    the multicomponent expression for the thermal conductivity.
%
help prandtl2

if nargin == 1 & isa(g,'solution')
   gas = g;
else
   gas = GRI30('Multi');
end

pr = zeros(31,31);
t = [];
xo2 = [];
io2 = speciesIndex(gas,'O2');
ih2 = speciesIndex(gas,'H2');

atm = oneatm;
t0 = cputime;
for i = 1:31
   t(i) = 300.0 + 100.0*i;
   for j = 1:31
      xo2(j) = 0.99*(j-1)/30.0;
      x = zeros(nSpecies(gas),1);      
      x(io2) = xo2(j);
      x(ih2) = 1.0 - xo2(j);
      set(gas,'T',t(i),'P',oneatm,'X',x);
      equilibrate(gas,'TP');
      pr(i,j) = viscosity(gas)*cp_mass(gas)/ ...
		thermalConductivity(gas);
      
   end
end
disp(['CPU time = ' num2str(cputime - t0)]);

% plot results

figure(1);
surf(xo2,t,pr);
xlabel('Elemental O/(O+H)');
ylabel('Temperature (K)');
zlabel('Prandtl Number');


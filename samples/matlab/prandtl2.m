function prandtl2(g)
% PRANDTL2  Prandtl number for an equilibrium H/O gas mixture.
%
%    This example does the same thing as prandtl1, but using
%    the multicomponent expression for the thermal conductivity.
%
% Keywords: transport, equilibrium, multicomponent transport, plotting

help prandtl2

if nargin == 1
   gas = g;
else
   gas = Solution('h2o2.yaml', 'ohmech', 'Multi');
end

pr = zeros(31,31);
xh2 = zeros(31,31);
visc = zeros(31,31);
lambda = zeros(31,31);
t = [];
xo2 = [];
io2 = speciesIndex(gas,'O2');
ih2 = speciesIndex(gas,'H2');

minT = minTemp(gas);
maxT = maxTemp(gas);
dT = (maxT - minT)/30.0;

t0 = cputime;
for i = 1:31
   t(i) = minT + dT*(i-1);
   for j = 1:31
      xo2(j) = 0.99*(j-1)/30.0;
      x = zeros(nSpecies(gas),1);
      x(io2) = xo2(j);
      x(ih2) = 1.0 - xo2(j);
      set(gas,'T',t(i),'P',oneatm,'X',x);
      equilibrate(gas,'TP');
      visc(i,j) = viscosity(gas);
      lambda(i,j) = thermalConductivity(gas);
      pr(i,j) = visc(i,j)*cp_mass(gas)/lambda(i,j);
      x = moleFractions(gas);
      xh2(i,j) = x(ih2);
   end
end
disp(['CPU time = ' num2str(cputime - t0)]);

% plot results

clf;
subplot(2,2,1);
surf(xo2,t,pr);
xlabel('Elemental O/(O+H)');
ylabel('Temperature (K)');
zlabel('Prandtl Number');

subplot(2,2,2);
surf(xo2,t,xh2);
xlabel('Elemental O/(O+H)');
ylabel('Temperature (K)');
zlabel('H_2 Mole Fraction');

subplot(2,2,3);
surf(xo2,t,visc);
xlabel('Elemental O/(O+H)');
ylabel('Temperature (K)');
zlabel('Viscosity');

subplot(2,2,4);
surf(xo2,t,lambda);
xlabel('Elemental O/(O+H)');
ylabel('Temperature (K)');
zlabel('Thermal Conductivity');

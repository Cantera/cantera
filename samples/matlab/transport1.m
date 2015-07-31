function transport1(g)
% TRANSPORT1  mixture-averaged transport properties
%
help transport1

if nargin == 1
   gas = g;
else
   gas = GRI30('Mix');
end

% set the state
%set(gas,'T',2300.0,'P',oneatm,'X','AR:1');

% = zeros(31,31);
pr = zeros(31,31);

t = [];
xo2 = [];
io2 = speciesIndex(gas,'O2');
ih2 = speciesIndex(gas,'H2');

atm = oneatm;
for i = 1:31
   t(i) = 300.0 + 100.0*i;
   for j = 1:31
      xo2(j) = (j-1)/30.0;
      x = zeros(nSpecies(gas),1);
      x(io2) = xo2(j);
      x(ih2) = 1.0 - xo2(j);
      set(gas,'T',t(i),'P',oneatm,'X',x);
      equilibrate(gas,'TP');
      pr(i,j) = viscosity(gas)*cp_mass(gas)/thermalConductivity(gas);
   end
end

% plot results

figure(1);
surf(xo2,t,pr);
xlabel('Elemental oxygen mole fraction');
ylabel('Temperature (K)');
zlabel('Prandtl Number');

function equil(g)
% EQUIL  a chemical equilibrium example.
%
%    This example computes the adiabatic flame temperature and
%    equilibrium composition for a methane/air mixture as a function of
%    equivalence ratio.
help equil

if nargin == 1
   gas = g;
else
   gas = Solution('gri30.yaml');
end

nsp = nSpecies(gas);

% find methane, nitrogen, and oxygen indices
ich4 = speciesIndex(gas,'CH4');
io2  = speciesIndex(gas,'O2');
in2  = speciesIndex(gas,'N2');

nPhis = 50;
phi = linspace(0.2, 2.70, nPhis);
tad(nPhis) = 0;
xeq(nsp,nPhis) = 0;

for i = 1:nPhis
   x = zeros(nsp,1);
   x(ich4,1) = phi(i);
   x(io2,1) = 2.0;
   x(in2,1) = 7.52;
   set(gas,'Temperature',300.0,'Pressure',101325.0,'MoleFractions',x);
   equilibrate(gas,'HP');
   tad(i) = temperature(gas);
   xeq(:,i) = moleFractions(gas);
end

% make plots
clf;
subplot(1,2,1);
plot(phi,tad);
xlabel('Equivalence Ratio');
ylabel('Temperature (K)');
title('Adiabatic Flame Temperature');

subplot(1,2,2);
semilogy(phi,xeq);
axis([phi(1) phi(50) 1.0e-14 1]);
%legend(speciesName(gas,1:nsp),1);
j = 10;
for k = 1:nsp
    text(phi(j),1.5*xeq(k,j),speciesName(gas,k))
    j = j + 2;
    if j > 46
        j = 10;
    end
end
xlabel('Equivalence Ratio');
ylabel('Mole Fraction');
title('Equilibrium Composition');

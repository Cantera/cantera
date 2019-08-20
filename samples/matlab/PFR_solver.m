function [F] = PFR_solver(x,init,gas,mdot,A_in,dAdx,k)

rho = init(1);
T = init(2);
Y = init(3:end);

if k==1
    A = A_in+k*dAdx*x;
elseif k==-1
    A = A_in+k*dAdx*x;
    dAdx = -dAdx;
else
    A = A_in+k*dAdx*x;
end

MW_mix = meanMolecularWeight(gas);
Ru = gasconstant();
R = Ru/MW_mix;
nsp = nSpecies(gas);
vx = mdot/(rho*A);
P = rho*R*T;

% the gas is set to the corresponding properties during each iteration of the ode loop
set(gas,'Temperature',T,'Density',rho,'MassFractions',Y);


MW = molecularWeights(gas);
h = enthalpies_RT(gas).*R.*T;
w = netProdRates(gas);
Cp = cp_mass(gas);
%--------------------------------------------------------------------------
%---F(1), F(2) and F(3:end) are the differential equations modelling the---
%---density, temperature and mass fractions variations along a plug flow---
%-------------------------reactor------------------------------------------
%--------------------------------------------------------------------------
F(1) = ((1-R/Cp)*((rho*vx)^2)*(1/A)*(dAdx)...
    + rho*R*sum(MW.*w.*(h-MW_mix*Cp*T./MW))/(vx*Cp) )...
    / (P*(1+vx^2/(Cp*T)) - rho*vx^2);

F(2) = (vx*vx/(rho*Cp))*F(1) + vx*vx*(1/A)*(dAdx)/Cp...
    - (1/(vx*rho*Cp))*sum(h.*w.*MW);

F(3:nsp+2) = w(1:nsp).*MW(1:nsp)./(rho*vx);

F = F';

end

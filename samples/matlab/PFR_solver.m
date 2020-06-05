function [F] = PFR_solver(x,soln_vector,gas,mdot,A_in,dAdx,k)

% This function defines the spatial derivatives for an ideal gas plug-flow 
% reactor, where the cross-sectional area and pressure are allowed to vary, 
% axially.  The model is set up by the example file 'Plug_Flow_Reactor.m', 
% which points the integrator to this function.  The integrator integrates the 
% derivatives spatially, to solve the density, temperature, and species mass 
% fraction profiles as a function of distance x.

rho = soln_vector(1);
T = soln_vector(2);
Y = soln_vector(3:end);

if k==1
    A = A_in+k*dAdx*x;
elseif k==-1
    A = A_in+k*dAdx*x;
    dAdx = -dAdx;
else
    A = A_in+k*dAdx*x;
end

% the gas is set to the corresponding properties during each iteration of the ode loop
set(gas,'Temperature',T,'Density',rho,'MassFractions',Y);

MW_mix = meanMolecularWeight(gas);
Ru = gasconstant();
R = Ru/MW_mix;
nsp = nSpecies(gas);
vx = mdot/(rho*A);
P = rho*R*T;

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

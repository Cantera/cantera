function isentropic(g)
% ISENTROPIC  isentropic, adiabatic flow example
%
%    In this example, the area ratio vs. Mach number curve is
%    computed for a hydrogen/nitrogen gas mixture.
%
help isentropic

if nargin == 1
   gas = g;
else
   gas = Solution('gri30.yaml');
end

% set the stagnation state
set(gas,'T',1200.0,'P',10.0*oneatm,'X','H2:1,N2:0.1');
s0 = entropy_mass(gas);
h0 = enthalpy_mass(gas);
p0 = pressure(gas);

mdot = 1;  % arbitrary

mach = [];
a = [];
i = 1;
amin = 1.e14;

% compute values for a range of pressure ratios
for r = 0.005:0.0025:0.995
   p = p0*r;

   % set the state using (p,s0)
   set(gas,'P',p,'S',s0);

   h = enthalpy_mass(gas);
   rho = density(gas);

   v2 = 2.0*(h0 - h);      %   h + V^2/2 = h0
   v = sqrt(v2);
   a(i) = mdot/(rho*v);    %   rho*v*A = constant

   if a(i) < amin
      amin = a(i);
   end
   mach(i) = v/soundspeed(gas);
   i = i + 1;
end

a = a/amin;

% plot results

clf;
plot(mach,a);
ylabel('Area Ratio');
xlabel('Mach Number');
title('Isentropic Flow: Area Ratio vs. Mach Number');

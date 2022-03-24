function f = flame(gas, left, flow, right)
%  FLAME - create a flame object.
%
%    gas     --  object representing the gas. This object will be used to
%                compute all required thermodynamic, kinetic, and transport
%                properties. The state of this object should be set
%                to an estimate of the gas state emerging from the
%                burner before calling StagnationFlame.
%
%    left    --  object representing the burner, which must be
%                created using function Inlet.
%
%    flow    --  object representing the flow, created with
%                function AxisymmetricFlow.
%
%    right   --  object representing the surface.
%

% Check input parameters
if nargin ~= 4
  error('wrong number of input arguments.');
end

if ~isIdealGas(gas)
  error('gas object must represent an ideal gas mixture.');
end
if ~isInlet(left)
  error('burner object of wrong type.');
end
if ~isFlow(flow)
  error('flow object of wrong type.');
end

flametype = 0;
if isSurface(right)
  flametype = 1;
elseif isInlet(right)
  flametype = 3;
end

% create the container object
f = Stack([left flow right]);

% set default initial profiles.
rho0 = density(gas);

% find the adiabatic flame temperature and corresponding
% equilibrium composition
equilibrate(gas, 'HP');
teq = temperature(gas);
yeq = massFractions(gas);

z1 = 0.2;
mdot0 = massFlux(left);
mdot1 = massFlux(right);
t0 = temperature(left);
if flametype == 0
  t1 = teq;
  mdot1 = -mdot0;
else
  t1 = temperature(right);
end
setProfile(f, 2, {'velocity', 'spread_rate'}, ...
           [0.0            1.0
            mdot0/rho0     -mdot1/rho0
            0.0            0.0]);
setProfile(f, 2, 'T', [0.0 z1 1.0; t0 2000.0 t1]);

for n = 1:nSpecies(gas)
  nm = speciesName(gas,n);
  if strcmp(nm,'H') | strcmp(nm,'OH') | strcmp(nm,'O') | strcmp(nm,'HO2')
    yint = 1.0*yeq(n);
  else
    yint = yeq(n);
  end
  if flametype == 3
    y1 = massFraction(right, n);
  else
    y1 = yeq(n);
  end
  setProfile(f, 2, nm, [0  z1   1
                        massFraction(left, n)   yint  y1]);
end

% set minimal grid refinement criteria
setRefineCriteria(f, 2, 10.0, 0.8, 0.8);

function flame = npflame_init(gas, left, flow, right, fuel, oxidizer, nuox)
%  FLAME - create a non-premixed flame object.
%
%    gas     --  object representing the gas. This object will be used to
%                compute all required thermodynamic, kinetic, and transport
%                properties. The state of this object should be set
%                to an estimate of the gas state emerging from the
%                burner before calling StagnationFlame.
%
%    left    --  object representing the left inlet, which must be
%                created using function Inlet.
%
%    flow    --  object representing the flow, created with
%                function AxisymmetricFlow.
%
%    right   --  object representing the right inlet, which must be
%                created using function Inlet.
%

% Check input parameters
if nargin ~= 7
    error('wrong number of input arguments.');
end

if ~isIdealGas(gas)
    error('gas object must represent an ideal gas mixture.');
end
if ~isInlet(left)
    error('left inlet object of wrong type.');
end
if ~isFlow(flow)
    error('flow object of wrong type.');
end
if ~isInlet(right)
    error('right inlet object of wrong type.');
end

% create the container object
flame = Stack([left flow right]);

% set default initial profiles.
rho0 = density(gas);

wt = molecularWeights(gas);

% find the fuel and oxidizer
ifuel = speciesIndex(gas,fuel);
ioxidizer = speciesIndex(gas,oxidizer);

s = nuox*wt(ioxidizer)/wt(ifuel);
y0f = massFraction(left,ifuel);
y0ox = massFraction(right,ioxidizer);
phi = s*y0f/y0ox;
zst = 1.0/(1.0 + phi);

% compute stoichiometric adiabatic flame temperature
nsp = nSpecies(gas);
tf = temperature(left);
tox = temperature(right);

yox = zeros(1, nsp);
yf = zeros(1, nsp);
ystoich = zeros(1, nsp);
for n = 1:nsp
    yox(n) = massFraction(right,n);
    yf(n) = massFraction(left,n);
    ystoich(n) = zst*yf(n) + (1.0 - zst)*yox(n);
end

set(gas,'T',temperature(left),'P',pressure(gas),'Y',ystoich);
equilibrate(gas,'HP');
teq = temperature(gas);
yeq = massFractions(gas);

% estimated strain rate
zz = z(flow);
dz = zz(end) - zz(1);
vleft = massFlux(left)/rho0;
vright = massFlux(right)/rho0;
a = (abs(vleft) + abs(vright))/dz;
diff = mixDiffCoeffs(gas);
f = sqrt(a/(2.0*diff(ioxidizer)));

x0 = massFlux(left)*dz/(massFlux(left) + massFlux(right));

nz = nPoints(flow);
zm = zeros(1,nz);
u = zeros(1,nz);
v = zeros(1,nz);
y = zeros(nz,nsp);
t = zeros(1,nz);
for j = 1:nz
    x = zz(j);
    zeta = f*(x - x0);
    zmix = 0.5*(1.0 - erf(zeta));
    zm(j) = zmix;
    u(j) = a*(x0 - zz(j));
    v(j) = a;
    if zmix > zst
        for n = 1:nsp
            y(j,n) = yeq(n) + (zmix - zst)*(yf(n) - yeq(n))/(1.0 - zst);
        end
        t(j) = teq + (tf - teq)*(zmix - zst)/(1.0 - zst);
    else
        for n = 1:nsp
            y(j,n) = yox(n) + zmix*(yeq(n) - yox(n))/zst;
        end
        t(j) = tox + zmix*(teq - tox)/zst;
    end
end

zrel = zz/dz;

setProfile(flame, 2, {'u', 'V'}, [zrel; u; v]);

setProfile(flame, 2, 'T', [zrel; t] );

for n = 1:nsp
    nm = speciesName(gas,n);
    setProfile(flame, 2, nm, [zrel; transpose(y(:,n))])
end

% set minimal grid refinement criteria
setRefineCriteria(flame, 2, 10.0, 0.99, 0.99);

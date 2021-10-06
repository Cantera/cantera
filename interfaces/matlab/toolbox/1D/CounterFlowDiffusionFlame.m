function flame = CounterFlowDiffusionFlame(left, flow, right, tp_f, tp_o, oxidizer)
% COUNTERFLOWDIFFUSIONFLAME  Create a counter flow diffusion flame stack.
% flame = CounterFlowDiffusionFlame(left, flow, right, tp_f, tp_o, oxidizer)
% :param left:
%     Object representing the left inlet, which must be
%     created using function :mat:func:`Inlet`.
% :param flow:
%     Object representing the flow, created with
%     function :mat:func:`AxisymmetricFlow`.
% :param right:
%     Object representing the right inlet, which must be
%     created using function :mat:func:`Inlet`.
% :param tp_f:
%     Object representing the fuel inlet gas, instance of class
%     :mat:func:`Solution`, and an ideal gas.
% :param tp_o:
%     Object representing the oxidizer inlet gas, instance of class
%     :mat:func:`Solution`, and an ideal gas.
% :param oxidizer:
%     String representing the oxidizer species. Most commonly O2.
% :return:
%     Instance of :mat:func:`Stack` object representing the left
%     inlet, flow, and right inlet.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input parameters
%

if nargin ~= 6
    error('CounterFlowDiffusionFlame expects six input arguments.');
end
if ~isIdealGas(tp_f)
    error('Fuel gas object must represent an ideal gas mixture.');
end
if ~isIdealGas(tp_o)
    error('Oxidizer gas object must represent an ideal gas mixture.');
end
if ~isInlet(left)
    error('Left inlet object of wrong type.');
end
if ~isFlow(flow)
    error('Flow object of wrong type.');
end
if ~isInlet(right)
    error('Right inlet object of wrong type.');
end
if ~ischar(oxidizer)
    error('Oxidizer name must be of format character.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the density of both fuel and oxidizer streams. To be used in
% determining velocity of each stream. Also get the temperature of both
% inlet streams.
%

rhof = density(tp_f);
rho0 = density(tp_o);
tf = temperature(left);
tox = temperature(right);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the species index of the oxidizer. To be used in determining initial
% strain rate.
%

ioxidizer = speciesIndex(tp_o, oxidizer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the stoichiometric mixture fraction. Needed for determining
% location of flame edges and composition. elMoles function used to
% calculate the number of moles of C, H, and O atoms in the fuel and
% oxidizer streams: elMoles = elementalMassFraction/element atomic weight.
% From this, the stoichiometric Air/Fuel ratio can be determined.
% 1 Mole of O needs 2 Moles of C and 0.5 Moles of H for stoichiometric
% conditions. The stoichiometric mixture fraction, Zst, is then calculated.
%

sFuel = elMoles(tp_f,'O')- 2*elMoles(tp_f,'C')- 0.5*elMoles(tp_f,'H');
sOx = elMoles(tp_o,'O')- 2*elMoles(tp_o,'C')- 0.5*elMoles(tp_o,'H');
phi = sFuel/sOx;
zst = 1.0/(1.0 - phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the stoichiometric mass fractions of each species. Use this to
% set the fuel gas object and calculate adiabatic flame temperature and
% equilibrium composition.
%

spec = speciesNames(tp_f); % Get all of the species names in gas object.
nsp = nSpecies(tp_f); % Get total number of species in gas object.
% Get the current mass fractions of both fuel and inlet streams.
yox = massFractions(tp_o);
yf = massFractions(tp_f);
ystoich_double = zeros(1, nsp); % Create empty vector for stoich mass frac.

for n = 1:nsp
    % Calculate stoichiometric mass fractions.
    ystoich_double(n) = zst*yf(n) + (1.0 - zst)*yox(n);
    % Convert mass fraction vector to string vector.
    ystoich_str{n} = num2str(ystoich_double(n));
    % Convert string vector to cell with SPECIES:MASS FRACTION format.
    y_stoich{n} = [spec{n}, ':', ystoich_str{n}];
end
% Initialize stoichiometric mass fraction cell with first SP:Y value.
ystoich = [y_stoich{1}];
for i = 2:nsp
    % Update cell to have format similar to N2:Yst,O2:Yst,...
    ystoich = [ystoich ',', y_stoich{i}];
end
% Set the fuel gas object as stoichiometric values and use equilibrate
% function to determine stoichiometric equilibrium temperature and mass
% fractions.
set(tp_f, 'T', tf, 'P', pressure(tp_f), 'Y', ystoich);
equilibrate(tp_f, 'HP');
teq = temperature(tp_f);
yeq = massFractions(tp_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the strain rate based on the inlet stream velocities and
% determine initial "guess" for mixture fraction based on mass flux ratio.
%

zz = gridPoints(flow);
dz = zz(end) - zz(1);
uleft = massFlux(left)/rhof;
uright = massFlux(right)/rho0;
a = (abs(uleft) + abs(uright))/dz;
diff = mixDiffCoeffs(tp_f);
f = sqrt(a/(2.0*diff(ioxidizer)));
x0num = sqrt(uleft*massFlux(left))*dz;
x0den = sqrt(uleft*massFlux(left)) + sqrt(uright*massFlux(right));
x0 = x0num/x0den;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial values of temperature and mass fraction of species in
% flame at each gridpoint. These values to be used for energy equation
% solution. Method is based on the Burke-Schumann model.
%

nz = nPoints(flow);
zm = zeros(1, nz);
u = zeros(1, nz);
v = zeros(1, nz);
y = zeros(nz, nsp);
t = zeros(1, nz);
for j = 1:nz
    x = zz(j);
    zeta = f*(x - x0);
    zmix = 0.5*(1.0 - erf(zeta)); % Mixture fraction in flame.
    zm(j) = zmix;
    u(j) = a*(x0 - zz(j)); % Axial velocity.
    v(j) = a; % Radial velocity.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the flame stack with the fuel inlet, flow object, and oxidizer
% inlet. Set the profile of the flame with the estimated axial velocities,
% radial velocities, temperature, and mass fractions calculated above.
%

flame = Stack([left flow right]);
setProfile(flame, 2, {'velocity', 'spread_rate'}, [zrel; u; v]);
setProfile(flame, 2, 'T', [zrel; t] );
for n = 1:nsp
    nm = speciesName(tp_f, n);
    setProfile(flame, 2, nm, [zrel; transpose(y(:,n))])
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define elMoles function
%

function moles = elMoles(tp, element)
% ELMOLES  Determine the elemental moles in a gas object per unit mass.
% moles = Moles(tp, element)
% :param tp:
%     Object representing the gas, instance of class :mat:func:`Solution`,
%     and an ideal gas. The state of this object should be set to an
%     estimate of the gas state before calling Moles.
% :param element:
%     String representing the element name.
% :return:
%     Elemental moles within a gas object per unit mass of mixture.
%     Units: kmol/kg
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input parameters
%

if nargin ~= 2
    error('elMoles expects two input arguments.');
end
if ~isIdealGas(tp)
    error('Gas object must represent an ideal gas mixture.');
end
if ~ischar(element)
    error('Element name must be of format character.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the moles per mass of mixture of an element within a gas
% object. The equation used is: elmoles = elMassFrac/Mel where elMassFrac
% is the elemental mass fraction within the gas object using the
% elementalMassFraction function; Mel is the atomic mass of the element.
%

elMassFrac = elementalMassFraction(tp, element);
eli = elementIndex(tp, element);
M = atomicMasses(tp);
Mel = M(eli);
moles = elMassFrac/Mel;
end

function flame = CounterFlowDiffusionFlame(left, flow, right, tp_f, tp_o, oxidizer)
    % Create a counter flow diffusion flame stack.
    % flame = CounterFlowDiffusionFlame(left, flow, right, tp_f, tp_o, oxidizer)
    %
    % :param left:
    %     Object representing the left inlet, which must be
    %     created using function :mat:class:`Inlet`.
    % :param flow:
    %     Object representing the flow, created with
    %     function :mat:class:`AxisymmetricFlow`.
    % :param right:
    %     Object representing the right inlet, which must be
    %     created using function :mat:class:`Inlet`.
    % :param tp_f:
    %     Object representing the fuel inlet gas, instance of class
    %     :mat:class:`Solution`, and an ideal gas.
    % :param tp_o:
    %     Object representing the oxidizer inlet gas, instance of class
    %     :mat:class:`Solution`, and an ideal gas.
    % :param oxidizer:
    %     String representing the oxidizer species. Most commonly O2.
    % :return:
    %     Instance of :mat:class:`Sim1D` object representing the left
    %     inlet, flow, and right inlet.
    %

    %% Check input parameters

    if nargin ~= 6
        error('CounterFlowDiffusionFlame expects six input arguments.');
    end

    if ~tp_f.isIdealGas
        error('Fuel gas object must represent an ideal gas mixture.');
    end

    if ~tp_o.isIdealGas
        error('Oxidizer gas object must represent an ideal gas mixture.');
    end

    if ~left.isInlet
        error('Left inlet object of wrong type.');
    end

    if ~flow.isFlow
        error('Flow object of wrong type.');
    end

    if ~right.isInlet
        error('Right inlet object of wrong type.');
    end

    if ~ischar(oxidizer)
        error('Oxidizer name must be of format character.');
    end

    %% Get the density of both fuel and oxidizer streams.
    % To be used in determining velocity of each stream. Also get the
    % temperature of both inlet streams.

    rhof = tp_f.D;
    rho0 = tp_o.D;
    tf = left.T;
    tox = right.T;

    %% Find the species index of the oxidizer.
    % To be used in determining initial strain rate.

    ioxidizer = tp_o.speciesIndex(oxidizer);

    %% Calculate the stoichiometric mixture fraction.
    % Needed for determining location of flame edges and composition.
    % elMoles function used to calculate the number of moles of C, H, and O
    % atoms in the fuel and oxidizer streams:
    %
    %    elMoles = elementalMassFraction/element atomic weight.
    %
    % From this, the stoichiometric Air/Fuel ratio can be determined.
    % 1 Mole of O needs 2 Moles of C and 0.5 Moles of H for stoichiometric
    % conditions. The stoichiometric mixture fraction, Zst, is then
    % calculated.

    sFuel = elMoles(tp_f, 'O') - 2 * elMoles(tp_f, 'C') - 0.5 * elMoles(tp_f, 'H');
    sOx = elMoles(tp_o, 'O') - 2 * elMoles(tp_o, 'C') - 0.5 * elMoles(tp_o, 'H');
    phi = sFuel / sOx;
    zst = 1.0 / (1.0 - phi);

    %% Compute the stoichiometric mass fractions of each species.
    % Use this to set the fuel gas object and calculate adiabatic flame
    % temperature and equilibrium composition.

    spec = tp_f.speciesNames; % Get all of the species names in gas object.
    nsp = tp_f.nSpecies; % Get total number of species in gas object.
    % Get the current mass fractions of both fuel and inlet streams.
    yox = tp_o.Y;
    yf = tp_f.Y;
    ystoich_double = zeros(1, nsp); % Create empty vector for stoich mass frac.

    for n = 1:nsp
        % Calculate stoichiometric mass fractions.
        ystoich_double(n) = zst * yf(n) + (1.0 - zst) * yox(n);
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
    tp_f.TPY = {tf, tp_f.P, ystoich};
    tp_f.equilibrate('HP');
    teq = tp_f.T;
    yeq = tp_f.Y;

    %% Estimate the strain rate.
    % Based on the inlet stream velocities and determine initial 'guess'
    % for mixture fraction based on mass flux ratio.

    zz = flow.gridPoints;
    dz = zz(end) - zz(1);
    mdotl = left.massFlux;
    mdotr = right.massFlux;
    uleft = mdotl / rhof;
    uright = mdotr / rho0;
    a = (abs(uleft) + abs(uright)) / dz;
    diff = tp_f.mixDiffCoeffs;
    f = sqrt(a / (2.0 * diff(ioxidizer)));
    x0num = sqrt(uleft * mdotl) * dz;
    x0den = sqrt(uleft * mdotr) + sqrt(uright * mdotr);
    x0 = x0num / x0den;

    %% Calculate initial values of temperature and mass fractions.
    % These values to be used for energy equation solution. Method is based
    % on the Burke-Schumann model.

    nz = flow.nPoints;
    zm = zeros(1, nz);
    u = zeros(1, nz);
    v = zeros(1, nz);
    y = zeros(nz, nsp);
    t = zeros(1, nz);

    for j = 1:nz
        x = zz(j);
        zeta = f * (x - x0);
        zmix = 0.5 * (1.0 - erf(zeta)); % Mixture fraction in flame.
        zm(j) = zmix;
        u(j) = a * (x0 - zz(j)); % Axial velocity.
        v(j) = a; % Radial velocity.

        if zmix > zst

            for n = 1:nsp
                y(j, n) = yeq(n) + (zmix - zst) * (yf(n) - yeq(n)) / (1.0 - zst);
            end

            t(j) = teq + (tf - teq) * (zmix - zst) / (1.0 - zst);
        else

            for n = 1:nsp
                y(j, n) = yox(n) + zmix * (yeq(n) - yox(n)) / zst;
            end

            t(j) = tox + zmix * (teq - tox) / zst;
        end

    end

    zrel = zz / dz;

    %% Create the flame stack.
    % Set the profile of the flame with the estimated axial velocities,
    % radial velocities, temperature, and mass fractions calculated above.
    flame = Sim1D([left flow right]);
    flame.setProfile(2, {'velocity', 'spread_rate'}, [zrel; u; v]);
    flame.setProfile(2, 'T', [zrel; t]);

    for n = 1:nsp
        nm = tp_f.speciesName(n);
        flame.setProfile(2, nm, [zrel; transpose(y(:, n))])
    end

end

%% Define elMoles function

function moles = elMoles(tp, element)
    % Determine the elemental moles in a gas object per unit mass.

    % Check input parameters
    if nargin ~= 2
        error('elMoles expects two input arguments.');
    end

    if ~tp.isIdealGas
        error('Gas object must represent an ideal gas mixture.');
    end

    if ~ischar(element)
        error('Element name must be of format character.');
    end

    % Calculate the moles per mass of mixture of an element within a gas
    % object. The equation used is:
    %
    %    elMoles = elMassFrac/Mel
    %
    % where elMassFrac is the elemental mass fraction within the gas object
    % using the elementalMassFraction function; Mel is the atomic mass of
    % the element.

    elMassFrac = tp.elementalMassFraction(element);
    eli = tp.elementIndex(element);
    M = tp.atomicMasses;
    Mel = M(eli);
    moles = elMassFrac / Mel;
end

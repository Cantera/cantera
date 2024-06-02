classdef Kinetics < handle
    % Kinetics Class ::
    %
    %     >> k = Kinetics(varargin)
    %
    % Class :mat:class:`Kinetics` represents kinetics managers, which manage
    % reaction mechanisms. The reaction mechanism attributes are specified in a
    % YAML file.
    %
    % Instances of class :mat:class:`Kinetics` are responsible for evaluating
    % reaction rates of progress, species production rates, and other
    % quantities pertaining to a reaction mechanism.
    %
    % :param varargin:
    %     Variable number of inputs consisting of the following:
    %       - ph:
    %           An instance of class :mat:class:`ThermoPhase` representing the
    %           phase in which reactions occur.
    %       - src:
    %           Input string of YAML file name when creating from file OR
    %           "clib" when called by the class constructors of :mat:class:`Solution`
    %           or :mat:class:`Interface`.
    %     Optional:
    %       - id:
    %           ID of the phase to import as specified in the input file.
    %       - neighbor1:
    %           Instance of class :mat:class:`ThermoPhase` or
    %           :mat:class:`Solution` representing the 1st neighboring phase.
    %       - neighbor2:
    %           Instance of class :mat:class:`ThermoPhase` or
    %           :mat:class:`Solution` representing the 2nd neighboring phase.
    %       - neighbor3:
    %           Instance of class :mat:class:`ThermoPhase` or
    %           :mat:class:`Solution` representing the 3rd neighboring phase.
    %       - neighbor4:
    %           Instance of class :mat:class:`ThermoPhase` or
    %           :mat:class:`Solution` representing the 4th neighboring phase.
    % :return:
    %      Instance of class :mat:class:`Kinetics`.

    properties (SetAccess = immutable)
        kinID % ID of the Kinetics object.
    end

    properties (SetAccess = protected)

        nPhases % Number of phases.

        nReactions % Number of reactions.

        nTotalSpecies % The total number of species.

        creationRates % Chemical reaction rates. Unit: kmol/m^3-s.

        destructionRates % Chemical destruction rates. Unit: kmol/m^3-s.

        deltaEnthalpy % Enthalpy of reaction. Unit: J/kmol.

        deltaStandardEnthalpy % Standard state enthalpy of reaction. Unit: J/kmol.

        deltaEntropy % Entropy of reaction. Unit: J/kmol-K.

        deltaStandardEntropy % Standard state entropy of reaction. Unit: J/kmol-K.

        deltaGibbs % Gibbs free energy of reaction. Unit: J/kmol-K.

        % Standard state gibbs free energy of reaction. Unit: J/kmol-K.
        deltaStandardGibbs

        % Equilibrium constants for all reactions. ::
        %
        %     >> k = kin.equilibriumConstants
        %
        % :return:
        %    A column vector of the equilibrium constants in concentration units for all
        %    reactions, calculated from the standard Gibbs free energy of reaction.
        equilibriumConstants

        forwardRateConstants % Forward reaction rate constants for all reactions.

        reverseRateConstants % Reverse reaction rate constants for all reactions.

        % Get the mass production rates of the species. ::
        %
        %     >> ydot = kin.massProdRate
        %
        % Evaluates the source term :math:`\dot{\omega}_k M_k /\rho`
        %
        % :param a:
        %     Instance of class :mat:class:`Kinetics` (or another
        %     object deriving from Kinetics)
        %     for which the ydots are desired.
        % :return:
        %     A vector of length nSpecies. Units: kg/s
        massProdRate

        netProdRates % Net chemical production rates for all species. Unit: kmol/m^3-s.

        % Forward rates of progress for all reactions. Unit: kmol/m^3-s.
        forwardRatesOfProgress

        % Reverse rates of progress for all reactions. Unit: kmol/m^3-s.
        reverseRatesOfProgress

        % Net rates of progress for all reactions. Unit: kmol/m^3-s.
        netRatesOfProgress

        reactionEquations % All reaction equations within the Kinetics object.

    end

    methods
        %% Kinetics Class Constructor

        function kin = Kinetics(varargin)

            ctIsLoaded;

            tmp = varargin{1};
            src = varargin{2};

            if ischar(tmp) & isnumeric(src)
                if strcmp(tmp, 'clib')
                    kin.kinID = ctFunc('soln_kinetics', src);
                    return
                end
            end

            ph = tmp.tpID;

            if nargin == 2
                id = '-';
            elseif nargin > 2
                id = varargin{3};
            end

            neighbours = {-1, -1, -1, -1};

            for i = 4:length(varargin)
                neighbours{i-3} = varargin{i}.tpID;
            end

            kin.kinID = ctFunc('kin_newFromFile', src, id, ph, neighbours{:});
        end

        %% Kinetics Class Destructor

        function delete(kin)
            % Delete the :mat:class:`Sim1D` object.

            if ~isa(kin, 'Solution') && ~isa(kin, 'Interface')
                ctFunc('kin_del', kin.kinID);
            end
        end

        %% Get scalar attributes

        function n = kineticsSpeciesIndex(kin, name, phase)
            % Get the species index of a species of a phase in the Kinetics class. ::
            %
            %     >> n = kin.kineticsSpeciesIndex(name, phase)
            %
            % :param name:
            %    String name or integer index of the species.
            % :param phase:
            %    String name or integer index of the phase.
            % :return:
            %    Index of the species.

            if nargin == 2
                phase = '';
            end
            n = ctFunc('kin_speciesIndex', kin.kinID, name, phase) + 1;
        end

        function n = multiplier(kin, irxn)
            % Get the multiplier for reaction rate of progress. ::
            %
            %     >> n = kin.multiplier(irxn)
            %
            % :param irxn:
            %    Integer reaction number for which the multiplier is desired.
            % :return:
            %    Multiplier of the rate of progress of reaction irxn.

            n = ctFunc('kin_multiplier', kin.kinID, irxn - 1);
        end

        function n = get.nPhases(kin)
            n = ctFunc('kin_nPhases', kin.kinID);
        end

        function n = get.nReactions(kin)
            n = ctFunc('kin_nReactions', kin.kinID);
        end

        function n = get.nTotalSpecies(kin)
            n = ctFunc('kin_nSpecies', kin.kinID);
        end

        function n = phaseIndex(kin, phase)
            % The index of a specific phase. ::
            %
            %     >> n = kin.phaseIndex(phase)
            %
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the phase.

            n = ctFunc('kin_phaseIndex', kin.kinID, phase) + 1;
        end

        function rxn = reactionEquation(kin, irxn)
            % Reaction equation of a reaction. ::
            %
            %   >> rxn = kin.reactionEqn(irxn)
            %
            % :param irxn:
            %    Integer index of the reaction.
            % :return:
            %    String reaction equation.

            rxn = ctString('kin_getReactionString', kin.kinID, irxn - 1);
        end

        %% Get reaction array attributes

        function n = reactantStoichCoeffs(kin, species, rxns)
            % Reactant stoichiometric coefficients. ::
            %
            %     >> n = kin.reactantStoichCoeffs(species, rxns)
            %
            % :param species:
            %    Species indices or names for which reactant stoichiometric
            %    coefficients should be retrieved.
            %    Optional argument; if specified, ``rxns`` must be specified as well.
            % :param rxns:
            %    Reaction indices for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, ``species`` must be specified as well.
            % :return:
            %    Returns a sparse matrix of all reactant stoichiometric
            %    coefficients if there are more than 1 species or reactions,
            %    otherwise returns an integer.
            %    The matrix elements ``nu(k, i)`` is the stoichiometric
            %    coefficient of species k as a reactant in reaction i.
            %    If ``species`` and ``rxns`` are specified, the matrix
            %    will contain only entries for the specified species
            %    and reactions. For example, ``kin.reactantStoichCoeffs(3,
            %    [1, 3, 5, 7])`` returns a sparse matrix containing only the
            %    coefficients for species 3 in reactions 1, 3, 5, and 7.

            nsp = kin.nTotalSpecies;
            nr = kin.nReactions;
            temp = sparse(nsp, nr);

            if nargin == 1
                krange = 1:nsp;
                irange = 1:nr;
            elseif nargin == 3
                if ischar(species)
                    krange = {species};
                else
                    krange = species;
                end
                irange = rxns;
            else
                error('stoichReactant requires 1 or 3 arguments.')
            end

            for k = krange

                for i = irange
                    if iscell(k)
                        kk = kin.kineticsSpeciesIndex(k{:});
                    elseif ischar(k)
                        kk = kin.kineticsSpeciesIndex(k);
                    else
                        kk = k;
                    end

                    t = ctFunc('kin_reactantStoichCoeff', kin.kinID, kk - 1, i - 1);

                    if t ~= 0.0
                        temp(kk, i) = t;
                    end

                end

            end

            if length(krange) > 1 || length(irange) > 1
                n = temp;
            else
                n = sum(full(temp), 'all');
            end
        end

        function n = productStoichCoeffs(kin, species, rxns)
            % Product stoichiometric coefficients. ::
            %
            %     >> n = kin.productStoichCoeffs(species, rxns)
            %
            % :param species:
            %    Species indices or names for which product stoichiometric
            %    coefficients should be retrieved.
            %    Optional argument; if specified, ``rxns`` must be specified as well.
            % :param rxns:
            %    Reaction indices for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, ``species`` must be specified as well.
            % :return:
            %    Returns a sparse matrix of all product stoichiometric
            %    coefficients if there are more than 1 species or reactions,
            %    otherwise returns an integer.
            %    The matrix elements ``nu(k, i)`` is the stoichiometric
            %    coefficient of species k as a product in reaction i.
            %    If ``species`` and ``rxns`` are specified, the matrix
            %    will contain only entries for the specified species
            %    and reactions. For example, ``kin.productStoichCoeffs(3,
            %    [1, 3, 5, 7])`` returns a sparse matrix containing only the
            %    coefficients for species 3 in reactions 1, 3, 5, and 7.

            nsp = kin.nTotalSpecies;
            nr = kin.nReactions;
            temp = sparse(nsp, nr);

            if nargin == 1
                krange = 1:nsp;
                irange = 1:nr;
            elseif nargin == 3
                if ischar(species)
                    krange = {species};
                else
                    krange = species;
                end
                irange = rxns;
            else
                error('stoichProduct requires 1 or 3 arguments.')
            end

            for k = krange

                for i = irange

                    if iscell(k)
                        kk = kin.kineticsSpeciesIndex(k{:});
                    else
                        kk = k;
                    end

                    t = ctFunc('kin_productStoichCoeff', kin.kinID, kk - 1, i - 1);

                    if t ~= 0.0
                        temp(kk, i) = t;
                    end

                end

            end

            if length(krange) > 1 || length(irange) > 1
                n = temp;
            else
                n = sum(full(temp), 'all');
            end
        end

        function n = netStoichCoeffs(kin, species, rxns)
            % Net stoichiometric coefficients. ::
            %
            %     >> n = kin.netStoichCoeffs(species, rxns)
            %
            % :param species:
            %    Species indices or names for which net stoichiometric
            %    coefficients should be retrieved.
            %    Optional argument; if specified, ``rxns`` must be specified as well.
            % :param rxns:
            %    Reaction indices for which net stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, ``species`` must be specified as well.
            % :return:
            %    Returns a sparse matrix of all net stoichiometric
            %    coefficients if there are more than 1 species or reactions,
            %    otherwise returns an integer.
            %    The matrix elements ``nu(k, i)`` is the stoichiometric
            %    coefficient of species k as a product in reaction i.
            %    If ``species`` and ``rxns`` are specified, the matrix
            %    will contain only entries for the specified species
            %    and reactions. For example, ``kin.netStoichCoeffs(3,
            %    [1, 3, 5, 7])`` returns a sparse matrix containing only the
            %    coefficients for species 3 in reactions 1, 3, 5, and 7.

            if nargin == 1
                n = kin.productStoichCoeffs - kin.reactantStoichCoeffs;
            elseif nargin == 3
                n = kin.productStoichCoeffs(species, rxns) - ...
                    kin.reactantStoichCoeffs(species, rxns);
            else
                error('stoichNet requires 1 or 3 arguments.');
            end

        end

        function cdot = get.creationRates(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getCreationRates', kin.kinID, nsp, pt);
            cdot = pt.Value;
        end

        function ddot = get.destructionRates(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getDestructionRates', kin.kinID, nsp, pt);
            ddot = pt.Value;
        end

        function n = isReversible(kin, i)
            % An array of flags indicating reversibility of a reaction. ::
            %
            %     >> n = kin.isReversible(i)
            %
            % :param i:
            %    Integer reaction number.
            % :return:
            %    True if reaction number i is reversible. false if irreversible.

            n = boolean(ctFunc('kin_isReversible', kin.kinID, i - 1));
        end

        function wdot = get.netProdRates(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getNetProductionRates', kin.kinID, nsp, pt);
            wdot = pt.Value;
        end

        function q = get.forwardRatesOfProgress(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getFwdRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = get.reverseRatesOfProgress(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getRevRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = get.netRatesOfProgress(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getNetRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function rxn = get.reactionEquations(kin)
            m = kin.nReactions;
            rxns = cell(1, m);

            for i = 1:m
                rxn{i} = kin.reactionEquation(i);
            end

        end

        function enthalpy = get.deltaEnthalpy(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getDelta', kin.kinID, 0, nr, pt);
            enthalpy = pt.Value;
        end

        function enthalpy = get.deltaStandardEnthalpy(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getDelta', kin.kinID, 3, nr, pt);
            enthalpy = pt.Value;
        end

        function entropy = get.deltaEntropy(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getDelta', kin.kinID, 2, nr, pt);
            entropy = pt.Value;
        end

        function entropy = get.deltaStandardEntropy(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getDelta', kin.kinID, 5, nr, pt);
            entropy = pt.Value;
        end

        function gibbs = get.deltaGibbs(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getDelta', kin.kinID, 1, nr, pt);
            gibbs = pt.Value;
        end

        function gibbs = get.deltaStandardGibbs(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getDelta', kin.kinID, 4, nr, pt);
            gibbs = pt.Value;
        end

        function k = get.equilibriumConstants(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getEquilibriumConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.forwardRateConstants(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getFwdRateConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.reverseRateConstants(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getRevRateConstants', kin.kinID, 0, nr, pt);
            k = pt.Value;
        end

        function ydot = get.massProdRate(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('kin_getSourceTerms', kin.kinID, nsp, pt);
            ydot = pt.Value;
        end

        %% Kinetics Set Methods

        function kin = setMultiplier(kin, irxn, v)
            % Set the multiplier for the reaction rate of progress. ::
            %
            %     >> kin.setMultiplier(irxn, v)
            %
            % :param irxn:
            %    Integer vector reaction numbers for which the multiplier
            %    should be set. Optional.
            % :param v:
            %    Value by which the reaction rate of progress should be
            %    multiplied.

            if nargin == 2
                v = irxn;
                n = kin.nReactions;
                irxn = (1:n);
            elseif nargin == 3
                n = length(irxn);
            else
                error('setMultiplier requires 2 or 3 arguments.')
            end

            for i = 1:n
                ctFunc('kin_setMultiplier', kin.kinID, irxn(i) - 1, v);
            end

        end

        function advanceCoverages(kin, dt)
            % Advance the surface coverages forward in time. ::
            %
            %     >> kin.advanceCoverages(dt)
            %
            % :param dt:
            %    Time interval by which the coverages should be advanced.

            ctFunc('kin_advanceCoverages', kin.kinID, dt);
        end

    end

end

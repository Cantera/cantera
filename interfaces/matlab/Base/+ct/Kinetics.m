classdef (Abstract) Kinetics < handle
    % Kinetics Class ::
    %
    %     >> k = ct.Kinetics(id)
    %
    % Retrieve instance of class :mat:class:`ct.Kinetics` associated with a
    % :mat:class:`ct.Solution` object. The constructor is called whenever a new
    % :mat:class:`ct.Solution` is instantiated and should not be used directly.
    %
    % Class :mat:class:`ct.Kinetics` represents kinetics managers, which manage
    % reaction mechanisms. The reaction mechanism attributes are specified in a
    % YAML file.
    %
    % Instances of class :mat:class:`ct.Kinetics` are responsible for evaluating
    % reaction rates of progress, species production rates, and other
    % quantities pertaining to a reaction mechanism.
    %
    % :param id:
    %     Integer ID of the solution holding the :mat:class:`ct.Kinetics` object.

    properties (SetAccess = immutable)
        kinID % ID of the Kinetics object.
    end

    properties (SetAccess = protected)

        nPhases % Number of phases.

        nReactions % Number of reactions.

        nTotalSpecies % The total number of species.

        creationRates % Creation rates for each species [kmol/m³/s for bulk phases]

        % Destruction rates for each species [kmol/m³/s for bulk phases]
        destructionRates

        deltaEnthalpy % Enthalpy of reaction [J/kmol]

        deltaStandardEnthalpy % Standard state enthalpy of reaction [J/kmol]

        deltaEntropy % Entropy of reaction [J/kmol/K]

        deltaStandardEntropy % Standard state entropy of reaction [J/kmol/K]

        deltaGibbs % Gibbs free energy of reaction [J/kmol]

        % Standard state Gibbs free energy of reaction [J/kmol]
        deltaStandardGibbs

        % Equilibrium constants in concentration units for all reactions, calculated
        % from the standard Gibbs free energy of reaction.
        equilibriumConstants

        forwardRateConstants % Forward reaction rate constants for all reactions.

        reverseRateConstants % Reverse reaction rate constants for all reactions.

        % Net chemical production rates for all species [kmol/m³/s for bulk phases]
        netProdRates

        % Forward rates of progress for all reactions [kmol/m³/s for bulk phases]
        forwardRatesOfProgress

        % Reverse rates of progress for all reactions [kmol/m³/s for bulk phases]
        reverseRatesOfProgress

        % Net rates of progress for all reactions [kmol/m³/s for bulk phases]
        netRatesOfProgress

        reactionEquations % All reaction equations within the Kinetics object.

    end

    methods
        %% Kinetics Class Constructor

        function kin = Kinetics(id)
            arguments
                id (1,1) double {mustBeInteger}
            end

            kin.kinID = ct.impl.call('mSol_kinetics', id);
        end

        %% Get scalar attributes

        function n = kineticsSpeciesIndex(obj, name)
            % Get the species index of a species of a phase in the Kinetics class. ::
            %
            %     >> n = kin.kineticsSpeciesIndex(name)
            %
            % :param name:
            %    String name or integer index of the species.
            % :return:
            %    Index of the species.

            n = ct.impl.call('mKin_kineticsSpeciesIndex', obj.kinID, name) + 1;
        end

        function n = multiplier(obj, irxn)
            % Get the multiplier for reaction rate of progress. ::
            %
            %     >> n = kin.multiplier(irxn)
            %
            % :param irxn:
            %    Integer reaction number for which the multiplier is desired.
            % :return:
            %    Multiplier of the rate of progress of reaction irxn.

            n = ct.impl.call('mKin_multiplier', obj.kinID, irxn - 1);
        end

        function n = get.nPhases(obj)
            n = ct.impl.call('mKin_nPhases', obj.kinID);
        end

        function n = get.nReactions(obj)
            n = ct.impl.call('mKin_nReactions', obj.kinID);
        end

        function n = get.nTotalSpecies(obj)
            n = ct.impl.call('mKin_nTotalSpecies', obj.kinID);
        end

        function n = phaseIndex(obj, phase)
            % The index of a specific phase. ::
            %
            %     >> n = kin.phaseIndex(phase)
            %
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the phase.

            n = ct.impl.call('mKin_phaseIndex', obj.kinID, phase) + 1;
        end

        function rxn = reactionEquation(obj, irxn)
            % Reaction equation of a reaction. ::
            %
            %   >> rxn = kin.reactionEquation(irxn)
            %
            % :param irxn:
            %    Integer index of the reaction.
            % :return:
            %    String reaction equation.
            r = ct.impl.call('mKin_reaction', obj.kinID, irxn - 1);
            rxn = ct.impl.getString('mRxn_equation', r);
        end

        %% Get reaction array attributes

        function n = reactantStoichCoeffs(obj, species, rxns)
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

            nsp = obj.nTotalSpecies;
            nr = obj.nReactions;
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
                        kk = obj.kineticsSpeciesIndex(k{:});
                    elseif ischar(k)
                        kk = obj.kineticsSpeciesIndex(k);
                    else
                        kk = k;
                    end

                    t = ct.impl.call('mKin_reactantStoichCoeff', obj.kinID, kk - 1, i - 1);

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

        function n = productStoichCoeffs(obj, species, rxns)
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

            nsp = obj.nTotalSpecies;
            nr = obj.nReactions;
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
                        kk = obj.kineticsSpeciesIndex(k{:});
                    else
                        kk = k;
                    end

                    t = ct.impl.call('mKin_productStoichCoeff', obj.kinID, kk - 1, i - 1);

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

        function n = netStoichCoeffs(obj, species, rxns)
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
                n = obj.productStoichCoeffs - obj.reactantStoichCoeffs;
            elseif nargin == 3
                n = obj.productStoichCoeffs(species, rxns) - ...
                    obj.reactantStoichCoeffs(species, rxns);
            else
                error('stoichNet requires 1 or 3 arguments.');
            end

        end

        function cdot = get.creationRates(obj)
            nsp = obj.nTotalSpecies;
            cdot = ct.impl.getArray('mKin_getCreationRates', nsp, obj.kinID);
        end

        function ddot = get.destructionRates(obj)
            nsp = obj.nTotalSpecies;
            ddot = ct.impl.getArray('mKin_getDestructionRates', nsp, obj.kinID);
        end

        function n = isReversible(obj, i)
            % An array of flags indicating reversibility of a reaction. ::
            %
            %     >> n = kin.isReversible(i)
            %
            % :param i:
            %    Integer reaction number.
            % :return:
            %    True if reaction number i is reversible. false if irreversible.

            n = ct.impl.call('mKin_isReversible', obj.kinID, i - 1) == 1;
        end

        function wdot = get.netProdRates(obj)
            nsp = obj.nTotalSpecies;
            wdot = ct.impl.getArray('mKin_getNetProductionRates', nsp, obj.kinID);
        end

        function q = get.forwardRatesOfProgress(obj)
            nr = obj.nReactions;
            q = ct.impl.getArray('mKin_getFwdRatesOfProgress', nr, obj.kinID);
        end

        function q = get.reverseRatesOfProgress(obj)
            nr = obj.nReactions;
            q = ct.impl.getArray('mKin_getRevRatesOfProgress', nr, obj.kinID);
        end

        function q = get.netRatesOfProgress(obj)
            nr = obj.nReactions;
            q = ct.impl.getArray('mKin_getNetRatesOfProgress', nr, obj.kinID);
        end

        function rxns = get.reactionEquations(obj)
            m = obj.nReactions;
            rxns = cell(1, m);

            for i = 1:m
                rxns{i} = obj.reactionEquation(i);
            end

        end

        function enthalpy = get.deltaEnthalpy(obj)
            nr = obj.nReactions;
            enthalpy = ct.impl.getArray('mKin_getDeltaEnthalpy', nr, obj.kinID);
        end

        function enthalpy = get.deltaStandardEnthalpy(obj)
            nr = obj.nReactions;
            enthalpy = ct.impl.getArray('mKin_getDeltaSSEnthalpy', nr, obj.kinID);
        end

        function entropy = get.deltaEntropy(obj)
            nr = obj.nReactions;
            entropy = ct.impl.getArray('mKin_getDeltaEntropy', nr, obj.kinID);
        end

        function entropy = get.deltaStandardEntropy(obj)
            nr = obj.nReactions;
            entropy = ct.impl.getArray('mKin_getDeltaSSEntropy', nr, obj.kinID);
        end

        function gibbs = get.deltaGibbs(obj)
            nr = obj.nReactions;
            gibbs = ct.impl.getArray('mKin_getDeltaGibbs', nr, obj.kinID);
        end

        function gibbs = get.deltaStandardGibbs(obj)
            nr = obj.nReactions;
            gibbs = ct.impl.getArray('mKin_getDeltaSSGibbs', nr, obj.kinID);
        end

        function k = get.equilibriumConstants(obj)
            nr = obj.nReactions;
            k = ct.impl.getArray('mKin_getEquilibriumConstants', nr, obj.kinID);
        end

        function k = get.forwardRateConstants(obj)
            nr = obj.nReactions;
            k = ct.impl.getArray('mKin_getFwdRateConstants', nr, obj.kinID);
        end

        function k = get.reverseRateConstants(obj)
            nr = obj.nReactions;
            k = ct.impl.getArray('mKin_getRevRateConstants', nr, obj.kinID, 'extraArgs', 0);
        end

        %% Kinetics Set Methods

        function kin = setMultiplier(obj, irxn, v)
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
                n = obj.nReactions;
                irxn = (1:n);
            elseif nargin == 3
                n = length(irxn);
            else
                error('setMultiplier requires 2 or 3 arguments.')
            end

            for i = 1:n
                ct.impl.call('mKin_setMultiplier', obj.kinID, irxn(i) - 1, v);
            end

        end

        function advanceCoverages(obj, dt)
            % Advance the surface coverages forward in time. ::
            %
            %     >> kin.advanceCoverages(dt)
            %
            % :param dt:
            %    Time interval by which the coverages should be advanced.

            ct.impl.call('mKin_advanceCoverages', obj.kinID, dt);
        end

    end

end

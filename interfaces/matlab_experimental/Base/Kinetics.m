classdef Kinetics < handle

    properties (SetAccess = immutable)
        kinID % ID of the Kinetics object.
    end

    properties (SetAccess = protected)

        nPhases % Number of phases.

        nReactions % Number of reactions.

        nTotalSpecies % The total number of species.

        creationRates % Chemical reaction rates. Unit: kmol/m^3-s.

        destructionRates % Chemical destruction rates. Unit: kmol/m^3-s.

        dH % enthalpy of reaction. Unit: J/kmol.

        dH_standard % standard state enthalpy of reaction. Unit: J/kmol.

        dS % entropy of reaction. Unit: J/kmol-K.

        dS_standard % standard state entropy of reaction. Unit: J/kmol-K.

        dG % gibbs free energy of reaction. Unit: J/kmol-K.

        dG_standard % standard state gibbs free energy of reaction. Unit: J/kmol-K.

        % Equilibrium constants for all reactions.
        %
        % k = kin.equilibriumConstants
        %
        % :return:
        %    A column vector of the equilibrium constants for all
        %    reactions. The vector has an entry for every reaction,
        %    whether reversible or not, but non-zero values occur only
        %    for the reversible reactions.
        equilibriumConstants

        forwardRateConstants % Forward reaction rate constants for all reactions.

        reverseRateConstants % Rever reaction rate constants for all reactions.

        % Get the mass production rates of the species.
        %
        % ydot = kin.massProdRate
        %
        % Evaluates the source term :math:`\dot{\omega}_k M_k /\rho`
        %
        % :param a:
        %     Instance of class :mat:class:`Kinetics` (or another
        %     object deriving from Kinetics)
        %     for which the ydots are desired.
        % :return:
        %     Returns a vector of length nSpecies. Units: kg/s
        massProdRate % Mass production rate of all species. Unit: kg/s.

        netProdRates % Net chemical production rates for all species. Unit: kmol/m^3-s.

        ropForward % Forward rates of progress for all reactions. Unit: kmol/m^3-s.

        ropReverse % Reverse rates of progress for all reactions. Unit: kmol/m^3-s.

        ropNet % Net rates of progress for all reactions. Unit: kmol/m^3-s.

        reactionEqns % All reaction equations within the Kinetics object.

    end

    methods
        %% Kinetics Class Constructor

        function kin = Kinetics(ph, src, id, n1, n2, n3, n4)
            % Kinetics Class
            %
            % k = Kinetics(ph, neighbor1, neighbor2, neighbor3, neighbor4)
            %
            % Class Kinetics represents kinetics managers, which are classes that manage
            % reaction mechanisms.  The reaction mechanism attributes are specified in a YAML file.
            %
            % Instances of class :mat:class:`Kinetics` are responsible for evaluating reaction rates
            % of progress, species production rates, and other quantities pertaining to
            % a reaction mechanism.
            %
            % :param ph:
            %     An instance of class :mat:class:`ThermoPhase` representing the phase
            %     in which reactions occur
            % :param src:
            %     Input string of YAML file name.
            % :param id:
            %     ID of the phase to import as specified in the input file. (optional)
            % :param neighbor1:
            %     Instance of class :mat:class:`ThermoPhase` or :mat:class:`Solution` representing a
            %     neighboring phase.
            % :param neighbor2:
            %     Instance of class :mat:class:`ThermoPhase` or :mat:class:`Solution` representing a
            %     neighboring phase.
            % :param neighbor3:
            %     Instance of class :mat:class:`ThermoPhase` or :mat:class:`Solution` representing a
            %     neighboring phase.
            % :param neighbor4:
            %     Instance of class :mat:class:`ThermoPhase` or :mat:class:`Solution` representing a
            %     neighboring phase.
            % :return:
            %      Instance of class :mat:class:`Kinetics`
            %
            checklib;

            % indices for bulk phases in a heterogeneous mechanism
            inb1 = -1;
            inb2 = -1;
            inb3 = -1;
            inb4 = -1;

            if nargin == 2
                id = '-';
            end

            % get the integer indices used to find the stored objects
            % representing the phases participating in the mechanism
            iph = ph.tpID;

            if nargin > 6
                inb4 = n4.tpID;
            end

            if nargin > 5
                inb3 = n3.tpID;
            end

            if nargin > 4
                inb2 = n2.tpID;
            end

            if nargin > 3
                inb1 = n1.tpID;
            end

            kin.kinID = callct('kin_newFromFile', src, id, ...
                                iph, inb1, inb2, inb3, inb4);
        end

        %% Kinetics Class Destructor

        function delete(kin)
            % Delete the kernel object

            callct('kin_del', kin.kinID);
        end

        %% Get scalar attributes

        function n = kineticsSpeciesIndex(kin, name, phase)
            % Get the species index of a species of a phase in the Kinetics class.
            %
            % :param name:
            %    String name of the species.
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the species.

            n = callct('kin_speciesIndex', kin.kinID, name, phase);
        end

        function n = multiplier(kin, irxn)
            % Get the multiplier for reaction rate of progress.
            %
            % :param irxn:
            %    Integer reaction number for which the multiplier is
            %    desired.
            % :return:
            %    Multiplier of the rate of progress of reaction irxn.

            n = callct('kin_multiplier', kin.kinID, irxn - 1);
        end

        function n = get.nPhases(kin)
            n = callct('kin_nPhases', kin.kinID);
        end

        function n = get.nReactions(kin)
            n = callct('kin_nReactions', kin.kinID);
        end

        function n = get.nTotalSpecies(kin)
            n = callct('kin_nSpecies', kin.kinID);
        end

        function n = phaseIndex(kin, phase)
            % The index of a specific phase.
            %
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the phase.

            n = callct('kin_phaseIndex', kin.kinID, phase);
        end

        function n = stoichReactant(kin, species, rxns)
            % Reactant stoichiometric coefficients.
            %
            % :param species:
            %    Species indices for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, ``rxns`` must be specified as well.
            % :param rxns:
            %    Reaction indicies for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, ``species`` must be specified as well.
            % :return:
            %    Returns a sparse matrix of all reactant stoichiometric
            %    coefficients. The matrix elements ``nu(k, i)`` is the
            %    stoichiometric coefficient of species k as a reactant in
            %    reaction i. If ``species`` and ``rxns`` are specified, the
            %    matrix will contain only entries for the specified species
            %    and reactions. For example, ``stoich_p(a, 3, [1, 3, 5,
            %    7])`` returns a sparse matrix containing only the
            %    coefficients for specis 3 in reactions 1, 3, 5, and 7.

            nsp = kin.nTotalSpecies;
            nr = kin.nReactions;
            temp = sparse(nsp, nr);

            if nargin == 1
                krange = 1:nsp;
                irange = 1:nr;
            elseif nargin == 3
                krange = species;
                irange = rxns;
            else error('stoichReactant requires 1 or 3 arguments.')
            end

            for k = krange

                for i = irange
                    t = callct('kin_reactantStoichCoeff', ...
                                kin.kinID, k - 1, i - 1);

                    if t ~= 0.0
                        temp(k, i) = t;
                    end

                end

            end

            n = temp;
        end

        function n = stoichProduct(kin, species, rxns)
            % Product stoichiometric coefficients.
            %
            % :param species:
            %    Species indices for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, ``rxns`` must be specified as well.
            % :param rxns:
            %    Reaction indicies for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, ``species`` must be specified as well.
            % :return:
            %    Returns a sparse matrix of all product stoichiometric
            %    coefficients.

            nsp = kin.nTotalSpecies;
            nr = kin.nReactions;
            temp = sparse(nsp, nr);

            if nargin == 1
                krange = 1:nsp;
                irange = 1:nr;
            elseif nargin == 3
                krange = species;
                irange = rxns;
            else error('stoichProduct requires 1 or 3 arguments.')
            end

            for k = krange

                for i = irange
                    t = callct('kin_productStoichCoeff', ...
                                kin.kinID, k - 1, i - 1);

                    if t ~= 0.0
                        temp(k, i) = t;
                    end

                end

            end

            n = temp;
        end

        function n = stoichNet(kin, species, rxns)
            % Net stoichiometric coefficients.
            %
            % :param species:
            %    Species indices for which net stoichiometric coefficients
            %    should be retrieved. Optional argument; if specified,
            %    "rxns" must be specified as well.
            % :param rxns:
            %    Reaction indicies for which net stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % :return:
            %    A sparse matrix of all net stoichiometric coefficients.

            if nargin == 1
                n = kin.stoichProduct - kin.stoichReactant;
            elseif nargin == 3
                n = kin.stoichProduct(species, rxns) - ...
                    kin.stoichReactant(species, rxns);
            else error('stoichNet requires 1 or 3 arguments.');
            end

        end

        %% Get reaction array attributes

        function cdot = get.creationRates(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getCreationRates', kin.kinID, nsp, pt);
            cdot = pt.Value;
        end

        function ddot = get.destructionRates(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDestructionRates', kin.kinID, nsp, pt);
            ddot = pt.Value;
        end

        function n = isReversible(kin, i)
            % An array of flags indicating reversibility of a reaction.
            %
            % n = kin.isReversible(i)
            %
            % :param i:
            %    Integer reaction number.
            % :return:
            %    1 if reaction number i is reversible. 0 if irreversible.

            n = callct('kin_isReversible', kin.kinID, i);
        end

        function wdot = get.netProdRates(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getNetProductionRates', kin.kinID, nsp, pt);
            wdot = pt.Value;
        end

        function q = get.ropForward(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getFwdRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = get.ropReverse(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getRevRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = get.ropNet(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getNetRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function rxn = reactionEqn(kin, irxn)
            % Reaction equation of a reaction
            %
            % rxn = kin.reactionEqn(irxn)
            %
            % :param irxn:
            %    Integer index of the reaction.
            % :return:
            %    String reaction equation.

            rxn = callct2('kin_getReactionString', kin.kinID, irxn - 1);
        end

        function rxn = get.reactionEqns(kin)
            m = kin.nReactions;
            rxns = cell(1, m);

            for i = 1:m
                rxn{i} = kin.reactionEqn(i);
            end

        end

        function enthalpy = get.dH(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 0, nr, pt);
            enthalpy = pt.Value;
        end

        function enthalpy = get.dH_standard(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 3, nr, pt);
            enthalpy = pt.Value;
        end

        function entropy = get.dS(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 2, nr, pt);
            entropy = pt.Value;
        end

        function entropy = get.dS_standard(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 5, nr, pt);
            entropy = pt.Value;
        end

        function gibbs = get.dG(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 1, nr, pt);
            gibbs = pt.Value;
        end

        function gibbs = get.dG_standard(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 4, nr, pt);
            gibbs = pt.Value;
        end

        function k = get.equilibriumConstants(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getEquilibriumConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.forwardRateConstants(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getFwdRateConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.reverseRateConstants(kin)
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getRevRateConstants', kin.kinID, 1, nr, pt);
            k = pt.Value;
        end

        function ydot = get.massProdRate(kin)
            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getSourceTerms', kin.kinID, nsp, pt);
            ydot = pt.Value;
        end

        %% Kinetics Set Methods

        function kin = setMultiplier(kin, irxn, v)
            % Set the multiplier for the reaction rate of progress.
            %
            % kin.setMultiplier(irxn, v)
            %
            % :param irxn:
            %    Integer vector reaction numbers for which the multiplier
            %    should be set. Optional.
            % :param v:
            %    Value by which the reaction rate of progress should be
            %    multiplied.
            %

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
                callct('kin_setMultiplier', kin.kinID, irxn(i) - 1, v);
            end

        end

        function advanceCoverages(kin, dt)
            % Advance the surface coveages forward in time
            %
            % kin.advanceCoverages(dt)
            %
            % :param dt:
            %    Time interval by which the coverages should be advanced.
            %

            callct('kin_advanceCoverages', kin.kinID, dt);
        end

    end

end

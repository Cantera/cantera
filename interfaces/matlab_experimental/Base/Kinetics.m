classdef Kinetics < handle
    % Kinetics Class.
    %
    % k = Kinetics(ph, neighbor1, neighbor2, neighbor3, neighbor4)
    %
    % Class Kinetics represents kinetics managers, which are classes that manage
    % reaction mechanisms.  The reaction mechanism attributes are specified in a YAML file.
    %
    % Instances of class :mat:func:`Kinetics` are responsible for evaluating reaction rates
    % of progress, species production rates, and other quantities pertaining to
    % a reaction mechanism.
    %
    % :param ph:
    %     An instance of class :mat:func:`ThermoPhase` representing the phase
    %     in which reactions occur
    % :param src:
    %     Input string of YAML file name.
    % :param id:
    %     ID of the phase to import as specified in the input file. (optional)
    % :param neighbor1:
    %     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
    %     neighboring phase.
    % :param neighbor2:
    %     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
    %     neighboring phase.
    % :param neighbor3:
    %     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
    %     neighboring phase.
    % :param neighbor4:
    %     Instance of class :mat:func:`ThermoPhase` or :mat:func:`Solution` representing a
    %     neighboring phase.
    % :return:
    %      Instance of class :mat:func:`Kinetics`
    %

    properties
        %% Scalar Attributes

        kinID % ID of the Kinetics object.
        kineticsSpeciesIndex
            % Get the species index of a species of a phase in the Kinetics class.
            %
            % :param name:
            %    String name of the species.
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the species.
            %
        multiplier
            % Get the multiplier for reaction rate of progress.
            %
            % :parameter irxn:
            %    Integer reaction number for which the multiplier is
            %    desired.
            % :return:
            %    Multiplier of the rate of progress of reaction irxn.
            %
        nPhases % Number of phases.
        nReactions % Number of reactions.
        nTotalSpecies % The total number of species.
        phaseIndex
            % The index of a specific phase.
            %
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the phase.
            %
        stoichReactant
            % Reactant stoichiometric coefficients.
            %
            % :parameter species:
            %    Species indices for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "rxns" must be specified as well.
            % :parameter rxns:
            %    Reaction indicies for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % :return:
            %    Returns a sparse matrix of all reactant stoichiometric
            %    coefficients. The matrix elements "nu(k, i)" is the
            %    stoichiometric coefficient of species k as a reactant in
            %    reaction i. If "species" and "rxns" are specified, the
            %    matrix will contain only entries for the specified species
            %    and reactions. For example, "stoich_p(a, 3, [1, 3, 5,
            %    7])" :returns a sparse matrix containing only the
            %    coefficients for specis 3 in reactions 1, 3, 5, and 7.
            %
        stoichProduct
            % Product stoichiometric coefficients.
            %
            % :parameter species:
            %    Species indices for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "rxns" must be specified as well.
            % :parameter rxns:
            %    Reaction indicies for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % :return:
            %    Returns a sparse matrix of all product stoichiometric
            %    coefficients.
            %
        stoichNet
            % Net stoichiometric coefficients.
            %
            % :parameter species:
            %    Species indices for which net stoichiometric coefficients
            %    should be retrieved. Optional argument; if specified,
            %    "rxns" must be specified as well.
            % :parameter rxns:
            %    Reaction indicies for which net stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % :return:
            %    A sparse matrix of all net stoichiometric coefficients.
            %

        %% Array Attributes
        creationRates % Chemical reaction rates. Unit: kmol/m^3-s.
        destructionRates % Chemical destruction rates. Unit: kmol/m^3-s.
        dH % enthalpy of reaction. Unit: J/kmol.
        dH_standard % standard state enthalpy of reaction. Unit: J/kmol.
        dS % entropy of reaction. Unit: J/kmol-K
        dS_standard % standard state entropy of reaction. Unit: J/kmol-K.
        dG % gibbs free energy of reaction. Unit: J/kmol-K.
        dG_standard % standard state gibbs free energy of reaction. Unit: J/kmol-K.
        equilibriumConstants
            % Equilibrium constants for all reactions.
            %
            % k = kin.equilibriumConstants
            %
            % :return:
            %    A column vector of the equilibrium constants for all
            %    reactions. The vector has an entry for every reaction,
            %    whether reversible or not, but non-zero values occur only
            %    for the reversible reactions.
            %
        forwardRateConstants % Forward reaction rate constants for all reactions.
        reverseRateConstants % Rever reaction rate constants for all reactions.
        isReversible
            % An array of flags indicating reversibility of a reaction.
            %
            % n = kin.isReversible(i)
            %
            % :parameter i:
            %    Integer reaction number.
            % :return:
            %    1 if reaction number i is reversible. 0 if irreversible.
            %
        massProdRate % Mass production rate of all species. Unit: kg/s.
        netProdRates % Net chemical production rates for all species. Unit: kmol/m^3-s.
        ropForward % Forward rates of progress for all reactions. Unit: kmol/m^3-s.
        ropReverse % Reverse rates of progress for all reactions. Unit: kmol/m^3-s.
        ropNet % Net rates of progress for all reactions. Unit: kmol/m^3-s.
        reactionEqn
            % Reaction equation of a reaction
            %
            % rxn = kin.reactionEqn(irxn)
            %
            % :parameter irxn:
            %    Integer index of the reaction.
            % :return:
            %    String reaction equation.
            %
        reactionEqns % All reaction equations within the Kinetics object.
    end

    methods
        %% Kinetics Class Constructor

        function kin = Kinetics(ph, src, id, n1, n2, n3, n4)

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
            if nargin > 3
                inb1 = n1.tpID;
                if nargin > 4
                    inb2 = n2.tpID;
                    if nargin > 5
                        inb3 = n3.tpID;
                        if nargin > 6
                            inb4 = n4.tpID;
                        end
                    end
                end
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

        function n = get.kineticsSpeciesIndex(kin, name, phase)
            % KINETICSSPECIESINDEX
            %
            % Get the species index in the Kinetics class.
            %
            % :param name:
            %    String name of the species.
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the species.
            %

            n = callct('kin_speciesIndex', kin.kinID, name, phase);
        end

        function n = get.multiplier(kin, irxn)
            % Get the multiplier for reaction rate of progress.
            %
            % :parameter irxn:
            %    Integer reaction number for which the multiplier is
            %    desired.
            % :return:
            %    Multiplier of the rate of progress of reaction irxn.

            n = callct('kin_multiplier', kin.kinID, irxn-1);
        end

        function n = get.nPhases(kin)
            % Get the number of phases.
            %
            % :return:
            %    Integer number of phases.

            n = callct('kin_nPhases', kin.kinID);
        end

        function n = get.nReactions(kin)
            % Get the number of reactions.
            %
            % :return:
            %    Integer number of reactions.

            n = callct('kin_nReactions', kin.kinID);
        end

        function n = get.nTotalSpecies(kin)
            % Get the total number of species.
            %
            % :return:
            %    Integer total number of species.

            n = callct('kin_nSpecies', kin.kinID);
        end

        function n = get.phaseIndex(kin, phase)
            % PHASEINDEX
            %
            % Get the index of a phase.
            %
            % :param phase:
            %    String name of the phase.
            % :return:
            %    Index of the phase.

            n = callct('kin_phaseIndex', kin.kinID, phase);
        end

        function n = get.stoichReactant(kin, species, rxns)
            % Get the reactant stoichiometric coefficients.
            %
            % :parameter species:
            %    Species indices for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "rxns" must be specified as well.
            % :parameter rxns:
            %    Reaction indicies for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % :return:
            %    Returns a sparse matrix of all reactant stoichiometric
            %    coefficients. The matrix elements "nu(k, i)" is the
            %    stoichiometric coefficient of species k as a reactant in
            %    reaction i. If "species" and "rxns" are specified, the
            %    matrix will contain only entries for the specified species
            %    and reactions. For example, "stoich_p(a, 3, [1, 3, 5,
            %    7])" :returns a sparse matrix containing only the
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
                                kin.kinID, k-1, i-1);
                    if t ~= 0.0
                        temp(k, i) = t;
                    end
                end
            end

            n = temp;
        end

        function n = get.stoichProduct(kin, species, rxns)
            % Get the product stoichiometric coefficients.
            %
            % :parameter species:
            %    Species indices for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "rxns" must be specified as well.
            % :parameter rxns:
            %    Reaction indicies for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % :return:
            %    Returns a sparse matrix of all product stoichiometric
            %    coefficients. The matrix elements "nu(k, i)" is the
            %    stoichiometric coefficient of species k as a product in
            %    reaction i. If "species" and "rxns" are specified, the
            %    matrix will contain only entries for the specified species
            %    and reactions. For example, "stoich_p(a, 3, [1, 3, 5,
            %    7])" :returns a sparse matrix containing only the
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
            else error('stoichProduct requires 1 or 3 arguments.')
            end

            for k = krange
                for i = irange
                    t = callct('kin_productStoichCoeff', ...
                                kin.kinID, k-1, i-1);
                    if t ~= 0.0
                        temp(k, i) = t;
                    end
                end
            end

            n = temp;
        end

        function n = get.stoichNet(kin, species, rxns)
            % Get the net stoichiometric coefficients.
            %
            % :parameter species:
            %    Species indices for which net stoichiometric coefficients
            %    should be retrieved. Optional argument; if specified,
            %    "rxns" must be specified as well.
            % :parameter rxns:
            %    Reaction indicies for which net stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % :return:
            %    A sparse matrix of all net stoichiometric coefficients.
            %    The matrix elements "nu(k, i)" is the stoichiometric
            %    coefficient of species k as a net in reaction.
            %    If "species" and "rxns" are specified, the matrix will
            %    contain only entries for the specified species and reactions.
            %    For example, "stoich_net(a, 3, [1, 3, 5, 7])" returns a
            %    sparse matrix containing only the coefficients for
            %    specis 3 in reactions 1, 3, 5, and 7.

            if nargin == 1
                n = kin.stoichProduct - kin.stoichReactant;
            elseif nargin == 3
                n = kin.stoichProduct(species, rxns) - kin.stoichReactant(species, rxns);
            else error('stoichNet requires 1 or 3 arguments.');
            end
        end

        %% Get reaction array attributes

        function cdot = get.creationRates(kin)
            % Get the chemical reaction rates.
            %
            % :return:
            %    A vector of the creation rates of all species. Unit: kmol/m^3-s.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getCreationRates', kin.kinID, nsp, pt);
            cdot = pt.Value;
        end

        function ddot = get.destructionRates(kin)
            % Get the chemical destruction rates.
            %
            % :return:
            %    A vector of the destruction rates of all species. Unit: kmol/m^3-s.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDestructionRates', kin.kinID, nsp, pt);
            ddot = pt.Value;
        end

        function n = get.isReversible(kin, i)
            % Get an array of flags indicating reversibility of a reaction.
            %
            % :parameter i:
            %    Integer reaction number.
            % :return:
            %    1 if reaction number i is reversible. 0 if irreversible.

            n = callct('kin_isReversible', kin.kinID, i);
        end

        function wdot = get.netProdRates(kin)
            % Get the net chemical production rates for all species.
            %
            % :return:
            %    A vector of the net production (creation-destruction)
            %    rates of all species. Unit: kmol/m^3-s.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getNetProductionRates', kin.kinID, nsp, pt);
            wdot = pt.Value;
        end

        function q = get.ropForward(kin)
            % Get the forward rates of progress for all reactions.
            %
            % :return:
            %    A column vector of the forward rates of progress for all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getFwdRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = get.ropReverse(kin)
            % Get the reverse rates of progress for all reactions.
            %
            % :return:
            %    A column vector of the reverse rates of progress for all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getRevRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = get.ropNet(kin)
            % Get the net rates of progress for all reactions.
            %
            % :return:
            %    A column vector of the net rates of progress for all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getNetRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function rxn = get.reactionEqn(kin, irxn)
            % Get the reaction equation of a reaction
            %
            % :parameter irxn:
            %    Integer index of the reaction.
            % :return:
            %    String reaction equation.

            output = callct2('kin_getReactionString', kin.kinID, irxn-1);
        end

        function rxn = get.reactionEqns(kin)
            % Get all reaction equations within the Kinetics class
            %
            % :return:
            %    Cell arrray of strings of the reaction equations.

            m = kin.nReactions;
            irxn = (1:m);
            rxns = cell(1, m);
            for i = 1:m
                rxns{i} = kin.reactionEqn(irxn(i)-1);
            end
        end

        function enthalpy = get.dH(kin)
            % Get the enthalpy of reaction for each reaction.
            %
            % :return:
            %    A vector of the enthalpy of reaction for each reaction.
            %    Unit: J/kmol.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 0, nr, pt);
            enthalpy = pt.Value;
        end

        function enthalpy = get.dH_standard(kin)
            % Get the standard state enthalpy of reaction for each reaction.
            %
            % :return:
            %    A vector of the standard state enthalpy of reaction for each reaction.
            %    Unit: J/kmol.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 3, nr, pt);
            enthalpy = pt.Value;
        end

        function entropy = get.dS(kin)
            % Get the entropy of reaction for each reaction.
            %
            % :return:
            %    A vector of the entropy of reaction for each reaction.
            %    Unit: J/kmol-K.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 2, nr, pt);
            entropy = pt.Value;
        end

        function entropy = get.dS_standard(kin)
            % Get the standard state entropy of reaction for each reaction.
            %
            % :return:
            %    A vector of the standard state entropy of reaction for each reaction.
            %    Unit: J/kmol-K.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 5, nr, pt);
            entropy = pt.Value;
        end

        function gibbs = get.dG(kin)
            % Get the Gibbs free energy of reaction for each reaction.
            %
            % :return:
            %    A vector of the Gibbs free energy of reaction for each
            %    reaction. Unit: J/kmol.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 1, nr, pt);
            gibbs = pt.Value;
        end

        function gibbs = get.dG_standard(kin)
            % Get the standard state Gibbs free energy of reaction for each reaction.
            %
            % :return:
            %    A vector of the standard state Gibbs free energy of reaction for each
            %    reaction. Unit: J/kmol.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getDelta', kin.kinID, 4, nr, pt);
            gibbs = pt.Value;
        end

        function k = get.equilibriumConstants(kin)
            % Get the equilibrium constants for all reactions.
            %
            % :return:
            %    A column vector of the equilibrium constants for all
            %    reactions. The vector has an entry for every reaction,
            %    whether reversible or not, but non-zero values occur only
            %    for the reversible reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getEquilibriumConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.forwardRateConstants(kin)
            % Get the forward reaction rate constants for all reactions.
            %
            % :return:
            %    A column vector of the forward rates constants of all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getFwdRateConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.reverseRateConstants(kin)
            % Get the reverse reaction rate constants for all reactions.
            %
            % :return:
            %    A column vector of the reverse rates constants of all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            callct('kin_getRevRateConstants', kin.kinID, 1, nr, pt);
            k = pt.Value;
        end

        function ydot = get.massProdRate(kin)
            % Get the mass production rates of the species.
            %
            % :return:
            %    A vector of the mass production rates.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('kin_getSourceTerms', kin.kinID, nsp, pt);
            massProdRate = pt.Value;
        end

        %% Kinetics Set Methods

        function kin = setMultiplier(kin, irxn, v)
            % Set the multiplier for the reaction rate of progress.
            %
            % kin.setMultiplier(irxn, v)
            %
            % :parameter irxn:
            %    Integer vector reaction numbers for which the multiplier
            %    should be set. Optional.
            % :parameter v:
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
                calct('kin_setMultiplier', kin.kinID, irxn(i)-1, v);
            end
        end

        function advanceCoverages(kin, dt)
            % Advance the surface coveages forward in time
            %
            % kin.advanceCoverages(dt)
            %
            % :parameter dt:
            %    Time interval by which the coverages should be advanced.
            %

            callct('kin_advanceCoverages', kin.kinID, dt);
        end

    end
end

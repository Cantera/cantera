classdef Kinetics < handle

    properties
        kinID
        Kc % equilibrium constant
        Kf % forward reaction rate
        Kr % reverse reaction rate
        dH % enthalpy of reaction
        dHss % standard state enthalpy of reaction
        dS % entropy of reaction
        dSss % standard state entropy of reaction
        dG % gibbs free energy of reaction
        dGss % standard state gibbs free energy of reaction
    end

    methods
        %% Kinetics class constructor

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
            kin.kinID = calllib(ct, 'kin_newFromFile', src, id, ...
                                 iph, inb1, inb2, inb3, inb4);
        end

        %% Utility methods

        function kinClear(kin)
            % Delete the kernel object

            calllib(ct, 'kin_del', kin.kinID);
        end

        %% Get scalar attributes

        function n = multiplier(kin, irxn)
            % Get the multiplier for reaction rate of progress.
            %
            % :parameter irxn:
            %    Integer reaction number for which the multiplier is
            %    desired.
            % :return:
            %    Multiplier of the rate of progress of reaction irxn.

            n = calllib(ct, 'kin_multiplier', kin.kinID, irxn-1);
        end

        function n = nReactions(kin)
            % Get the number of reactions.
            %
            % :return:
            %    Integer number of reactions

            n = calllib(ct, 'kin_nReactions', kin.kinID);
        end

        function n = nTotalSpecies(kin)
            % Get the total number of species.
            %
            % :return:
            %    Integer total number of species.

            n = calllib(ct, 'kin_nSpecies', kin.kinID);
        end

        function n = stoichReactant(kin, species, rxns)
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
            else error('stoich_r requires 1 or 3 arguments.')
            end

            for k = krange
                for i = irange
                    t = calllib(ct, 'kin_reactantStoichCoeff', ...
                                kin.kinID, k-1, i-1);
                    if t ~= 0.0
                        temp(k, i) = t;
                    end
                end
            end

            n = temp;
        end

        function n = stoichProduct(kin, species, rxns)
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
            else error('stoich_p requires 1 or 3 arguments.')
            end

            for k = krange
                for i = irange
                    t = calllib(ct, 'kin_productStoichCoeff', ...
                                kin.kinID, k-1, i-1);
                    if t ~= 0.0
                        temp(k, i) = t;
                    end
                end
            end

            n = temp;
        end

        function n = stoichNet(kin, species, rxns)
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
                n = stoich_p(kin)-stoich_r(kin);
            elseif nargin == 3
                n = stoich_p(a, species, rxns) - stoich_r (a, species, rxns);
            else error('stoich_net requires 1 or 3 arguments.');
            end
        end

        %% Get reaction array attributes

        function cdot = creationRates(kin)
            % Get the chemical reaction rates.
            %
            % :return:
            %    A vector of the creation rates of all species. Unit: kmol/m^3-s.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getCreationRates', kin.kinID, nsp, pt);
            cdot = pt.Value;
        end

        function ddot = destructionRates(kin)
            % Get the chemical destruction rates.
            %
            % :return:
            %    A vector of the destruction rates of all species. Unit: kmol/m^3-s.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDestructionRates', kin.kinID, nsp, pt);
            ddot = pt.Value;
        end

        function n = isReversible(kin, i)
            % Get an array of flags indicating reversibility of a reaction.
            %
            % :parameter i:
            %    Integer reaction number.
            % :return:
            %    1 if reaction number i is reversible. 0 if irreversible.

            n = calllib(ct, 'kin_isReversible', kin.kinID, i);
        end

        function wdot = netProdRates(kin)
            % Get the net chemical production rates for all species.
            %
            % :return:
            %    A vector of the net production (creation-destruction)
            %    rates of all species. Unit: kmol/m^3-s.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getNetProductionRates', kin.kinID, nsp, pt);
            wdot = pt.Value;
        end

        function q = ropForward(kin)
            % Get the forward rates of progress for all reactions.
            %
            % :return:
            %    A column vector of the forward rates of progress for all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getFwdRateOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = ropReverse(kin)
            % Get the reverse rates of progress for all reactions.
            %
            % :return:
            %    A column vector of the reverse rates of progress for all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getRevRateOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function q = ropNet(kin)
            % Get the net rates of progress for all reactions.
            %
            % :return:
            %    A column vector of the net rates of progress for all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getNetRatesOfProgress', kin.kinID, nr, pt);
            q = pt.Value;
        end

        function rxn = reactionEqn(kin, irxn)
            % Get the reaction equation of a reaction
            %
            % :parameter irxn:
            %    Optional. Integer or vector of reaction numbers.
            % :return:
            %    String or cell arrray of strings of the reaction
            %    equations.

            if nargin == 1
                m = kin.nReactions;
                n = 1
                irxn = (n:m)';
            elseif nargin == 2
                if isa(irxn, 'double')
                    [m, n] = size(irxn);
                else
                    error('reaction numbers must be numeric');
                end
            end

            rxn = cell(m, n);
            for i = 1:m
                for j = 1:n
                    buflen = calllib(ct, 'kin_getReactionString', kin.kinID, ...
                                     irxn - 1, 0, '');
                    if buflen > 0
                            aa = char(zeros(1, buflen));
                            [~, aa] = calllib(ct, 'kin_getReactionString', ...
                                              kin.kinID, irxn - 1, buflen, aa);
                            rxn{i, j} = aa;
                    end
                end
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
            calllib(ct, 'kin_getDelta', kin.kinID, 0, nr, pt);
            enthalpy = pt.Value;
        end

        function enthalpy = get.dHss(kin)
            % Get the standard state enthalpy of reaction for each reaction.
            %
            % :return:
            %    A vector of the standard state enthalpy of reaction for each reaction.
            %    Unit: J/kmol.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDelta', kin.kinID, 3, nr, pt);
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
            calllib(ct, 'kin_getDelta', kin.kinID, 2, nr, pt);
            entropy = pt.Value;
        end

        function entropy = get.dSss(kin)
            % Get the standard state entropy of reaction for each reaction.
            %
            % :return:
            %    A vector of the standard state entropy of reaction for each reaction.
            %    Unit: J/kmol-K.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDelta', kin.kinID, 5, nr, pt);
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
            calllib(ct, 'kin_getDelta', kin.kinID, 1, nr, pt);
            gibbs = pt.Value;
        end

        function gibbs = get.dGss(kin)
            % Get the standard state Gibbs free energy of reaction for each reaction.
            %
            % :return:
            %    A vector of the standard state Gibbs free energy of reaction for each
            %    reaction. Unit: J/kmol.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDelta', kin.kinID, 4, nr, pt);
            gibbs = pt.Value;
        end

        function k = get.Kc(kin)
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
            calllib(ct, 'kin_getEquilibriumConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.Kf(kin)
            % Get the forward reaction rate constants for all reactions.
            %
            % :return:
            %    A column vector of the forward rates constants of all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getFwdRateConstants', kin.kinID, nr, pt);
            k = pt.Value;
        end

        function k = get.Kr(kin)
            % Get the reverse reaction rate constants for all reactions.
            %
            % :return:
            %    A column vector of the reverse rates constants of all
            %    reactions.

            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getRevRateConstants', kin.kinID, 1, nr, pt);
            k = pt.Value;
        end

        function massProdRate = ydot(kin)
            % Get the mass production rates of the species.
            %
            % :return:
            %    A vector of the mass production rates.

            nsp = kin.nTotalSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getSourceTerms', kin.kinID, nsp, pt);
            massProdRate = pt.Value;
        end

        %% Set attributes

        function n = setMultiplier(kin, irxn, v)
            % Set the multiplier for the reaction rate of progress.
            %
            % :parameter irxn:
            %    Integer of vector reaction numbers for which the
            %    multiplier should be set. Optional.
            % :parameter v:
            %    Value by which the reaction rate of progress should be
            %    multiplied.

            if nargin == 2
                v = irxn;
                nr = kin.nReactions;
                irxn = (1:nr)';
                n = 1;
            elseif nargin == 3
                [nr, n] = size(irxn);
            else
                error('setMultiplier requires 2 or 3 arguments.')
            end

            for i = 1:nr
                for j = 1:n
                    calllib(ct, 'kin_setMultiplier', kin.kinID, ...
                            irxn(i, j)-1, v);
                end
            end
        end

        function advanceCoverages(kin, dt)
            % Advance the surface coveages forward in time
            %
            % :parameter dt:
            %    Time interval by which the coverages should be advanced.

            calllib(ct, 'kin_advanceCoverages', kin.kinID, dt);
        end

    end
end

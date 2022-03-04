classdef Kinetics < handle

    properties
        owner
        kin_id
        Kc % equilibrium constant
        Kf % forward reaction rate
        Kr % reverse reaction rate
        dH % enthalpy of reaction
        dS % entropy of reaction
        dG % gibbs free energy of reaction
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
            kin.owner = 1;
            % get the integer indices used to find the stored objects
            % representing the phases participating in the mechanism
            iph = ph.tp_id;
            if nargin > 3
                inb1 = n1.tp_id;
                if nargin > 4
                    inb2 = n2.tp_id;
                    if nargin > 5
                        inb3 = n3.tp_id;
                        if nargin > 6
                            inb4 = n4.tp_id;
                        end
                    end
                end
            end
            kin.kin_id = calllib(ct, 'kin_newFromFile', src, id, ...
                                 iph, inb1, inb2, inb3, inb4);
        end

        %% Utility methods

        function kin_clear(kin)
            % Delete the kernel object

            checklib;
            calllib(ct, 'kin_del', kin.kin_id);
        end

        %% Get scalar attributes

        function n = isReversible(kin, i)
            % Get an array of flags indicating reversibility of a reaction.
            %
            % parameter i:
            %    Integer reaction number.
            % return:
            %    1 if reaction number i is reversible. 0 if irreversible.

            checklib;
            n = calllib(ct, 'kin_isReversible', kin.kin_id, i);
        end

        function n = multiplier(kin, irxn)
            % Get the multiplier for reaction rate of progress.
            %
            % parameter irxn:
            %    Integer reaction number for which the multiplier is
            %    desired.
            % return:
            %    Multiplier of the rate of progress of reaction irxn.

            checklib;
            n = calllib(ct, 'kin_multiplier', kin.kin_id, irxn-1);
        end

        function n = nReactions(kin)
            % Get the number of reactions.
            %
            % return:
            %    Integer number of reactions

            checklib;
            n = calllib(ct, 'kin_nReactions', kin.kin_id);
        end

        function n = nSpecies2(kin)
            % Get the total number of species.
            %
            % return:
            %    Integer total number of species.

            checklib;
            n = calllib(ct, 'kin_nSpecies', kin.kin_id);
        end

        function n = stoich_r(kin, species, rxns)
            % Get the reactant stoichiometric coefficients.
            %
            % parameter species:
            %    Species indices for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "rxns" must be specified as well.
            % parameter rxns:
            %    Reaction indicies for which reactant stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % return:
            %    Returns a sparse matrix of all reactant stoichiometric
            %    coefficients. The matrix elements "nu(k, i)" is the
            %    stoichiometric coefficient of species k as a reactant in
            %    reaction i. If "species" and "rxns" are specified, the
            %    matrix will contain only entries for the specified species
            %    and reactions. For example, "stoich_p(a, 3, [1, 3, 5,
            %    7])" returns a sparse matrix containing only the
            %    coefficients for specis 3 in reactions 1, 3, 5, and 7.

            checklib;
            nsp = kin.nSpecies;
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
                                kin.kin_id, k-1, i-1);
                    if t ~= 0.0
                        temp(k, i) = t;
                    end
                end
            end

            n = temp;
        end

        function n = stoich_p(kin, species, rxns)
            % Get the product stoichiometric coefficients.
            %
            % parameter species:
            %    Species indices for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "rxns" must be specified as well.
            % parameter rxns:
            %    Reaction indicies for which product stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % return:
            %    Returns a sparse matrix of all product stoichiometric
            %    coefficients. The matrix elements "nu(k, i)" is the
            %    stoichiometric coefficient of species k as a product in
            %    reaction i. If "species" and "rxns" are specified, the
            %    matrix will contain only entries for the specified species
            %    and reactions. For example, "stoich_p(a, 3, [1, 3, 5,
            %    7])" returns a sparse matrix containing only the
            %    coefficients for specis 3 in reactions 1, 3, 5, and 7.

            checklib;
            nsp = kin.nSpecies;
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
                                kin.kin_id, k-1, i-1);
                    if t ~= 0.0
                        temp(k, i) = t;
                                            end
                end
            end

            n = temp;
        end

        function n = stoich_net(kin, species, rxns)
            % Get the net stoichiometric coefficients.
            %
            % parameter species:
            %    Species indices for which net stoichiometric coefficients
            %    should be retrieved. Optional argument; if specified,
            %    "rxns" must be specified as well.
            % parameter rxns:
            %    Reaction indicies for which net stoichiometric
            %    coefficients should be retrieved. Optional argument; if
            %    specified, "species" must be specified as well.
            % return:
            %    Returns a sparse matrix of all net stoichiometric
            %    coefficients. The matrix elements "nu(k, i)" is the
            %    stoichiometric coefficient of species k as a net in
            %    reaction i. If "species" and "rxns" are specified, the
            %    matrix will contain only entries for the specified species
            %    and reactions. For example, "stoich_net(a, 3, [1, 3, 5,
            %    7])" returns a sparse matrix containing only the
            %    coefficients for specis 3 in reactions 1, 3, 5, and 7.

            checklib;
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
            % return:
            %    Returns a vector of the creation rates of all species. If
            %    the output is not assigned to a variable, a bar graph is
            %    produced. Unit: kmol/m^3-s.

            checklib;
            nsp = kin.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getCreationRates', kin.kin_id, nsp, pt);
            cdot = pt.Value;
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Creation Rates')
%                 bar(q)
%                 xlabel('Species Number')
%                 ylabel('Creation Rate [kmol/m^3-s]')
%                 title('Species Chemical Reaction Rates')
%             end
        end

        function ddot = destructionRates(kin)
            % Get the chemical destruction rates.
            %
            % return:
            %    Returns a vector of the destruction rates of all species.
            %    If the output is not assigned to a variable, a bar graph
            %    is produced. Unit: kmol/m^3-s.

            checklib;
            nsp = kin.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDestructionRates', kin.kin_id, nsp, pt);
            ddot = pt.Value;
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Destruction Rates')
%                 bar(q)
%                 xlabel('Species Number')
%                 ylabel('Destruction Rate [kmol/m^3-s]')
%                 title('Species Chemical Reaction Rates')
%             end
        end

        function wdot = netProdRates(kin)
            % Get the net chemical production rates for all species.
            %
            % return:
            %    Returns a vector of the net production (creation-destruction)
            %    rates of all species. If the output is not assigned to a
            %    variable, a bar graph is produced. Unit: kmol/m^3-s.

            checklib;
            nsp = kin.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getNetProductionRates', kin.kin_id, nsp, pt);
            wdot = pt.Value;
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Production Rates')
%                 bar(q)
%                 xlabel('Species Number')
%                 ylabel('Net Production Rate [kmol/m^3-s]')
%                 title('Species Net Chemical Reaction Rates')
%             end
        end

        function q = rop_f(kin)
            % Get the forward rates of progress for all reactions.
            %
            % return:
            %    Returns a column vector of the forward rates of progress
            %    for all reactions. If this function is called without
            %    argument, a bar graph is produced.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getFwdRateOfProgress', kin.kin_id, nr, pt);
            q = pt.Value;
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Rates of Progress')
%                 bar(q)
%                 xlabel('Reaction Number')
%                 ylabel('Forward Rate of Progress [kmol/m^3]')
%                 title('Forward Rates of Progress')
%             end
        end

        function q = rop_r(kin)
            % Get the reverse rates of progress for all reactions.
            %
            % return:
            %    Returns a column vector of the reverse rates of progress
            %    for all reactions. If this function is called without
            %    argument, a bar graph is produced.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getRevRateOfProgress', kin.kin_id, nr, pt);
            q = pt.Value;
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Rates of Progress')
%                 bar(q)
%                 xlabel('Reaction Number')
%                 ylabel('Reverse Rate of Progress [kmol/m^3]')
%                 title('Reverse Rates of Progress')
%             end
        end

        function q = rop(kin)
            % Get the forward and reverse rates of progress.
            %
            % return:
            %    Returns an I x 2 array of reaction rates of progress,
            %    where I is the number of reactions. The first column
            %    contains the forward rates of progress, and the second
            %    column the reverse rates. If this function is called
            %    without arguments, a bar graph is produced.

            checklib;
            f = rop_f(kin);
            r = rop_r(kin);
            q = [f, r];
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Rates of Progress')
%                 bar(q)
%                 xlabel('Reaction Number')
%                 ylabel('Rate of Progress [kmol/m^3]')
%                 title('Rates of Progress')
%                 legend('Forward', 'Reverse')
%             end
        end

        function q = rop_net(kin)
            % Get the net rates of progress for all reactions.
            %
            % return:
            %    Returns a column vector of the net rates of progress
            %    for all reactions. If this function is called without
            %    argument, a bar graph is produced.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getNetRatesOfProgress', kin.kin_id, nr, pt);
            q = pt.Value;
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Rates of Progress')
%                 bar(q)
%                 xlabel('Reaction Number')
%                 ylabel('Net Rate of Progress [kmol/m^3]')
%                 title('Net Rates of Progress')
%             end
        end

        function rxn = reactionEqn(kin, irxn)
            % Get the reaction equation of a reaction
            %
            % parameter irxn:
            %    Optional. Integer or vector of reaction numbers.
            % return:
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
                    buflen = calllib(ct, 'kin_getReactionString', kin.kin_id, ...
                                     irxn - 1, 0, '');
                    if buflen > 0
                            aa = char(zeros(1, buflen));
                            [~, aa] = calllib(ct, 'kin_getReactionString', ...
                                              kin.kin_id, irxn - 1, buflen, aa);
                            rxn{i, j} = aa;
                    end
                end
            end
        end

        function enthalpy = get.dH(kin)
            % Get the enthalpy of reaction for each reaction.
            %
            % return:
            %    Returns a vector of the enthalpy of reaction for each
            %    reaction. Unit: J/kmol.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDelta', kin.kin_id, 0, nr, pt);
            enthalpy = pt.Value;
        end

        function entropy = get.dS(kin)
            % Get the entropy of reaction for each reaction.
            %
            % return:
            %    Returns a vector of the entropy of reaction for each
            %    reaction. Unit: J/kmol-K.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDelta', kin.kin_id, 2, nr, pt);
            entropy = pt.Value;
        end

        function gibbs = get.dG(kin)
            % Get the Gibbs free energy of reaction for each reaction.
            %
            % return:
            %    Returns a vector of the Gibbs free energy of reaction for
            %    each reaction. Unit: J/kmol.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getDelta', kin.kin_id, 1, nr, pt);
            gibbs = pt.Value;
        end

        function k = get.Kc(kin)
            % Get the equilibrium constants for all reactions.
            %
            % return:
            %    Returns a column vector of the equilibrium constants for
            %    all reactions. The vector has an entry for every reaction,
            %    whether reversible or not, but non-zero values occur only
            %    for the reversible reactions. If the output is not
            %    assigned to a variable, a bar graph is produced.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getEquilibriumConstants', kin.kin_id, nr, pt);
            k = pt.Value;
%             if nargout == 0
%                 figure
%                 set(gcf, 'Name', 'Equilibrium Constants')
%                 bar(k)
%                 xlabel('Reaction Number')
%                 ylabel('log_{10} Kc [kmol,m, s]')
%                 title('Equilibrium Constants')
%             end
        end

        function k = get.Kf(kin)
            % Get the forward reaction rate constants for all reactions.
            %
            % return:
            %    Returns a column vector of the forward rates constants of
            %    all of the reactions.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getFwdRateConstants', kin.kin_id, nr, pt);
            k = pt.Value;
        end

        function k = get.Kr(kin)
            % Get the reverse reaction rate constants for all reactions.
            %
            % return:
            %    Returns a column vector of the reverse rates constants of
            %    all of the reactions.

            checklib;
            nr = kin.nReactions;
            xx = zeros(1, nr);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getRevRateConstants', kin.kin_id, 1, nr, pt);
            k = pt.Value;
        end

        function massProdRate = ydot(kin)
            % Get the mass production rates of the species.
            %
            % return:
            %    Returns a vector of the mass production rates.

            checklib;
            nsp = kin.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'kin_getSourceTerms', kin.kin_id, nsp, pt);
            massProdRate = pt.Value;
        end

        %% Set attributes

        function n = setMultiplier(kin, irxn, v)
            % Set the multiplier for the reaction rate of progress.
            %
            % parameter irxn:
            %    Integer of vector reaction numbers for which the
            %    multiplier should be set. Optional.
            % parameter v:
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
                    calllib(ct, 'kin_setMultiplier', kin.kin_id, ...
                            irxn(i, j)-1, v);
                end
            end
        end

        function advanceCoverages(kin, dt)
            % Advance the surface coveages forward in time
            %
            % parameter dt:
            %    Time interval by which the coverages should be advanced.

            calllib(ct, 'kin_advanceCoverages', kin.kin_id, dt);
        end

    end
end

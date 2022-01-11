classdef ThermoPhase < handle

    properties
        tp_owner
        tp_id
        T % temperature
        P % pressure
        D % density
        X % mole fractions
        Y % mass fractions
        H % enthalpy
        S % entropy
        U % internal energy
        G % Gibbs free energy
        V % specific volume
        basis
    end

    properties (Dependent)
        DP
        DPX
        DPY
        HP
        HPX
        HPY
        PV
        PVX
        PVY
        SH
        SHX
        SHY
        SP
        SPX
        SPY
        ST
        STX
        STY
        SV
        SVX
        SVY
        TD
        TDX
        TDY
        TH
        THX
        THY
        TP
        TPX
        TPY
        TV
        TVX
        TVY
        UV
        UVX
        UVY
        UP
        UPX
        UPY
        VH
        VHX
        VHY
    end

    methods
        %% ThermoPhase class constructor

        function tp = ThermoPhase(src, id)
            checklib;
            if nargin > 2
                error('ThermoPhase expects 1 or 2 input arguments.');
            end
            if nargin == 1
                id = '-';
            end
            tp.tp_owner = 1;
            tp.tp_id = calllib(ct, 'thermo_newFromFile', src, id);
            tp.basis = 'molar';
        end

        %% Utility methods

        function display(tp, threshold)
            % Display thermo properties

            if nargin < 2 || ~isnumeric(threshold)
                threshold = 1e-14;
            end
            calllib(ct, 'thermo_print', tp.tp_id, 1, threshold);
        end

        function tp_clear(tp)
            % Delete the kernel object.

            checklib;
            calllib(ct, 'thermo_del', tp.tp_id);
        end

        function tp = set.basis(tp, b)
            % Set basis of thermodynamic properties
            % Default is molar

            if strcmp(b, 'mole') || strcmp(b, 'molar') ...
               || strcmp(b, 'Mole') || strcmp(b, 'Molar')
                tp.basis = 'molar';
            elseif strcmp(b, 'mass') || strcmp(b, 'Mass')
                tp.basis = 'mass';
            else error("Basis must be mass or molar")
            end
        end

        %% PhaseGet single methods

        function amu = atomicMasses(tp)
            % Get the atomic masses of the elements.
            %
            % return:
            %    Vector of element atomic masses. Unit: kg/kmol

            checklib;
            nel = tp.nElements;
            aa = zeros(1, nel);
            pt = libpointer('doublePtr', aa);
            calllib(ct, 'thermo_getAtomicWeights', ...
                    tp.tp_id, nel, pt);
            amu = pt.Value;
        end

        function e = charges(tp)
            % Get the array of species charges.
            %
            % return:
            %    Vector of species charges. Unit: elem. charge

            checklib;
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(ct, 'thermo_getCharges', ...
                    tp.tp_id, nsp, pt);
            e = pt.Value;
        end

        function k = elementIndex(tp, name)
            % Get the index of an element given the name.
            % The index is an integer assigned to each element in sequence
            % as it is read in from the input file.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.
            %
            % parameter name:
            %    String or cell array of elements whose index is requested
            % return:
            %    Integer number of elements in the phase.

            checklib;
            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);
                for i = 1:m
                    for j = 1:n
                        k(i, j) = calllib(ct, 'thermo_elementIndex', ...
                                          tp.tp_id, name{i, j}) + 1;
                        if k(i, j) > 1e3
                            warning(['Element ', name{i, j}, ...
                                   ' does not exist in the phase']);
                            k(i, j) = -1;
                        end
                    end
                end
            elseif ischar(name)
                k = calllib(ct, 'thermo_elementIndex', ...
                            tp.tp_id, name) + 1;
                if k > 1e3
                    warning(['Element ', name, ...
                           ' does not exist in the phase']);
                    k = -1;
                end
            else
                error('name must be either a string or cell array of strings')
            end
        end

        function elMassFrac = elementalMassFraction(tp, element)
            % Determine the elemental mass fraction in gas object.
            checklib;

            % Check input parameters.
            if nargin ~= 2
                error('elementalMassFraction expects two input arguments.');
            end
            if ~isIdealGas(tp)
                error('Gas object must represent an ideal gas mixture.');
            end
            if ~ischar(element)
                error('Element name must be of format character.');
            end

            % Calculate the elemental mass fraction in a gas object using
            % the following equation:
            %
            %    elMassFrac = sum of nAtoms(k, m)*Mel(m)*Y(k)/mw(k)
            %
            % where nAtoms(k, m) is the number of atoms of element 'm', in
            % species 'k'; Mel(m) is the atomic weight of element 'm'; Y(k)
            % is the mass fraction of species 'k'; and mw(k) is the
            % molecular weight of species 'k'.

            n = tp.nSpecies;
            spec = tp.speciesNames;
            eli = tp.elementIndex(element);
            M = tp.atomicMasses;
            Mel = M(eli);
            MW = tp.MolecularWeights;
            % Initialize the element mass fraction as zero.
            elMassFrac = 0.0;
            % Use loop to perform summation of elemental mass fraction over all species.
            for i = 1:n
                natoms(i) = tp.nAtoms(spec{i}, element);
                yy(i) = tp.massFraction(spec{i});
                elMassFrac = elMassFrac + (natoms(i)*Mel*yy(i))/MW(i);
            end
        end

        function mmw = meanMolecularWeight(tp)
            % Get the mean molecular weight.
            %
            % return:
            %    Double mean molecular weight. Unit: kg/kmol

            checklib;
            mmw = calllib(ct, 'thermo_meanMolecularWeight', tp.tp_id);
        end

        function density = molarDensity(tp)
            % Get the molar basis density in kmol/m^3.
            checklib;
            density = calllib(ct, 'thermo_molarDensity', tp.tp_id);
        end

        function mw = MolecularWeights(tp)
            % Get the array of molecular weights of all species.
            %
            % return:
            %   Vector of species molecular weights. Unit: kg/kmol

            checklib;
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(ct, 'thermo_getMolecularWeights', ...
                    tp.tp_id, nsp, pt);
            mw = pt.Value;
        end

        function n = nAtoms(tp, species, element)
            % Get the number of atoms of an element in a species.
            %
            % parameter k:
            %    String species name or integer species number.
            % parameter m:
            %    String element name or integer element number.
            % return:
            %    Integer number of atoms of the element in the species.
            if nargin == 3
                if ischar(species)
                    k = tp.speciesIndex(species);
                else k = species;
                end
                if k < 0
                    n = -1;
                    return
                end
                if ischar(element)
                    m = tp.elementIndex(element);
                else m = element;
                end
                if m < 0
                    n = -1;
                    return
                end
                n = calllib(ct, 'thermo_nAtoms', tp.tp_id, k-1, m-1);
            else
                error('Two input arguments required.')
            end
        end

        function nel = nElements(tp)
            % Get the number of elements.
            %
            % return:
            %    Integer number of elements in the phase.

            checklib;
            nel = calllib(ct, 'thermo_nElements', tp.tp_id);
        end

        function nsp = nSpecies(tp)
            % Get the number of species.
            %
            % return:
            %    Integer number of species in the phase.

            checklib;
            nsp = calllib(ct, 'thermo_nSpecies', tp.tp_id);
        end

        function k = speciesIndex(tp, name)
            % Get the index of a species given the name.
            % The index is an integer assigned to each species in sequence
            % as it is read in from the input file.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.
            %
            % parameter name:
            %    String or cell array of species whose index is requested.
            % return:
            %    Integer number of species in the phase.

            checklib;
            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);
                for i = 1:m
                    for j = 1:n
                        k(i, j) = calllib(ct, 'thermo_speciesIndex', ...
                                          tp.tp_id, name{i, j}) + 1;
                        if k(i, j) > 1e3
                            warning(['Species ', name{i, j}, ...
                                   ' does not exist in the phase']);
                            k(i, j) = -1;
                        end
                    end
                end
            elseif ischar(name)
                k = calllib(ct, 'thermo_speciesIndex', ...
                            tp.tp_id, name) + 1;
                if k > 1e3
                    warning(['Species ', name, ...
                           ' does not exist in the phase.']);
                    k = -1;
                end
            else
                error('Name must be either a string or cell array of strings.')
            end
        end

        function nm = speciesName(tp, k)
            % Get the name of a species given the index.
            %
            % parameter k:
            %    Scalar or array of integer species index.
            % return:
            %    Cell array of strings species name.
            [m, n] = size(k);
            nm = cell(m, n);
            for i = 1:m
                for j = 1:n
                    ksp = k(i, j) - 1;
                    buflen = calllib(ct, 'thermo_getSpeciesName', ...
                                     tp.tp_id, ksp, 0, '');
                    if buflen > 0
                        aa = char(zeros(1, buflen));
                        [~, aa] = calllib(ct, 'thermo_getSpeciesName', ...
                                                tp.tp_id, ksp, buflen, aa);
                        nm{i, j} = aa;
                    end
                end
            end
        end

        function n = speciesNames(tp)
            % Get all species names.
            n = tp.speciesName(1:tp.nSpecies);
        end

        function temperature = get.T(tp)
            % Get the temperature.
            %
            % return:
            %    Double temperature. Unit: K

            checklib;
            temperature = calllib(ct, 'thermo_temperature', tp.tp_id);
        end

        function pressure = get.P(tp)
            % Get the pressure.
            %
            % return:
            %    Double pressure. Unit: Pa

            checklib;
            pressure = calllib(ct, 'thermo_pressure', tp.tp_id);
        end

        function density = get.D(tp)
            % Get the mass basis density in kg/m^3.
            checklib;
            density = calllib(ct, 'thermo_density', tp.tp_id);
        end

        function volume = get.V(tp)
            % Get the specific volume depending on the basis.
            %
            % return:
            %    Density depending on the basis. Units:
            %    m^3/kmol (molar) m^3/kg (mass).
            volume = 1/tp.D;
        end

        function moleFractions = get.X(tp)
            % Get the mole fractions of all species.
            %
            % return:
            %    Vector of species mole fractions.

            checklib;
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'thermo_getMoleFractions', ...
                    tp.tp_id, nsp, pt);
            moleFractions = pt.Value;

            % if no output argument is specified, a bar plot is produced.
            if nargout == 0
                figure
                set(gcf, 'Name', 'Mole Fractions')
                bar(moleFractions)
                xlabel('Species Number')
                ylabel('Mole Fraction')
                title('Species Mole Fractions')
            end
        end

        function x = moleFraction(tp, species)
            % Get the mole fraction of one or a list of species.
            %
            % parameter species:
            %    String or cell array of species whose mole fraction is
            %    requested.
            % return:
            %    Scalar or vector of species mole fractions.

            xarray = tp.X;
            if isa(species, 'char')
                k = tp.speciesIndex(species);
                if  k > 0
                    x = xarray(k);
                else error("species not found.");
                end
            elseif isa(species, 'cell')
                n = length(species);
                x = zeros(1, n);
                for j = 1:n
                    k = tp.speciesIndex(species{j});
                    if k > 0
                        x(j) = xarray(k);
                    else error("species not found.");
                    end
                end
            end
        end

        function massFractions = get.Y(tp)
            % Get the mass fractions of all species.
            %
            % return:
            %   Vector of species mass fractions.

            checklib;
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(ct, 'thermo_getMassFractions', ...
                    tp.tp_id, nsp, pt);
            massFractions = pt.Value;

            % If no output argument is specified, a bar plot is produced.
            if nargout == 0
                figure
                set(gcf, 'Name', 'Mole Fractions')
                bar(massFractions)
                xlabel('Species Number')
                ylabel('Mole Fraction')
                title('Species Mole Fractions')
            end
        end

        function y = massFraction(tp, species)
            % Get the mass fraction of one or a list of species.
            %
            % parameter species:
            %    String or cell array of species whose mass fraction is
            %    requested.
            % return:
            %    Scalar or vector of species mass fractions.

            yy = tp.Y;
            if isa(species, 'char')
                k = tp.speciesIndex(species);
                if  k > 0
                    y = yy(k);
                else error("Error: species not found.")
                end
            elseif isa(species, 'cell')
                n = length(species);
                y = zeros(1, n);
                for j = 1:n
                    k = tp.speciesIndex(species{j});
                    if k > 0
                        y(j) = yy(k);
                    else error("Error: species not found.")
                    end
                end
            end
        end

        %% ThermoGet single methods

        function mu = chemical_potentials(tp)
            % Get the chemical potentials of the species.
            %
            % return:
            %    Vector of species chemical potentials. Unit: J/kmol.

            checklib;
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'thermo_chemPotentials', ...
                    tp.tp_id, nsp, pt);
            mu = pt.Value;
        end

        function c = cv(tp)
            % Get the specific heat at constant volume.
            %
            % return:
            %    Specific heat of the mixture at constant volume depending
            %    on the basis. Units: J/kmol-K (molar) J/kg-K (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                c = calllib(ct, 'thermo_cv_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                c = calllib(ct, 'thermo_cv_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function c = cp(tp)
            % Get the specific heat at constant pressure.
            %
            % return:
            %    Specific heat of the mixture at constant pressure depending
            %    on the basis. Units: J/kmol-K (molar) J/kg-K (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                c = calllib(ct, 'thermo_cp_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                c = calllib(ct, 'thermo_cp_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function d = critDensity(tp)
            % Get the critical density.
            %
            % return:
            %    Critical density. Unit: K.

            checklib;
            d = calllib(ct, 'thermo_critDensity', tp.tp_id);
        end

        function p = critPressure(tp)
            % Get the critical pressure.
            %
            % return:
            %    Critical temperature. Unit: Pa.

            checklib;
            p = calllib(ct, 'thermo_critPressure', tp.tp_id);
        end

        function t = critTemperature(tp)
            % Get the critical temperature.
            %
            % return:
            %    Critical temperature. Unit: K.
            checklib;
            t = calllib(ct, 'thermo_critTemperature', tp.tp_id);
        end

        function v = electricPotential(tp)
            % Get the electric potential
            %
            % return:
            %    Electric potential of the phase. Unit: V.

            checklib;
            v = calllib(ct, 'thermo_electricPotential', tp.tp_id);
        end

        function e = eosType(tp)
            % Get the type of the equation of state
%             checklib;
%             buflen = calllib(ct, 'thermo_getEosType', tp.tp_id, 0, '');
%             if buflen > 0
%                 aa = char(zeros(1, buflen));
%                 [~, aa] = calllib(ct, 'thermo_getEosType', ...
%                                         tp.tp_id, buflen, aa);
%             end
%             e = aa;
            e = 'IdealGas';
        end

        function v = isIdealGas(tp)
            % Get a flag indicating whether the phase is an ideal gas.

            if strcmp(tp.eosType, 'IdealGas')
                v = 1;
            else
                v = 0;
            end
            v = 1;
        end

        function b = isothermalCompressibility(tp)
            % Get the isothermal compressibility
            %
            % return:
            %    Isothermal compressibility. Unit: 1/Pa.

            checklib;
            b = calllib(ct, 'thermo_isothermalCompressibility', tp.tp_id);
        end

        function t = maxTemp(tp)
            % Get the maximum temperature of the parameter fits.
            %
            % return:
            %    Vector of maximum temperatures of all species.

            checklib;
            t = calllib(ct, 'thermo_maxTemp', tp.tp_id, -1);
        end

        function t = minTemp(tp)
            % Get the minimum temperature of the parameter fits.
            %
            % return:
            %    Vector of minimum temperatures of all species.

            checklib;
            t = calllib(ct, 'thermo_minTemp', tp.tp_id, -1);
        end

        function p = P_sat(tp, t)
            % Get the saturation pressure for a given temperature.
            %
            % parameter t:
            %    Temperature. Unit: K.
            % return:
            %    Saturation pressure for temperature t. Unit: Pa.

            checklib;
            p = calllib(ct, 'thermo_satPressure', tp.tp_id, t);
        end

        function p = refPressure(tp)
            % Get the reference pressure.
            %
            % return:
            %    Reference pressure. Unit: Pa.

            checklib;
            p = calllib(ct, 'thermo_refPressure', tp.tp_id, -1);
        end

        function c = soundspeed(tp)
            % Get the speed of sound
            % If the phase is an ideal gas, the speed of sound is
            % calculated by:
            %
            %    c = sqrt(gamma * R * T);
            %
            % where gamma is the ratio of specific heat, R is the specific
            % gas constant, and T is the temperature. If the phase is not
            % an ideal gas, the speed of sound is calculated by:
            %
            %    c = sqrt(
            %
            % return:
            %    The speed of sound. Unit: m/s

            checklib;
            if tp.isIdealGas
                tp.basis = 'mass';
                gamma = tp.cp/tp.cv;
                wtm = tp.meanMolecularWeight;
                r = 8314.4621/wtm;
                c = sqrt(gamma * r *tp.T);
            else
                rho0 = tp.D;
                p0 = tp.P;
                s0 = tp.S;
                rho1 = 1.001 * rho0;
                tp.DS = {rho1, s0};
                p1 = tp.P;
                dpdrho_s = (p1 - p0) / (rho1 - rho0);
                c = sqrt(dpdrho_s);
            end
        end

        function a = thermalExpansionCoeff(tp)
            % Get the thermal expansion coefficient.
            %
            % return:
            %    Thermal expansion coefficient. Unit: 1/K.

            checklib;
            a = calllib(ct, 'thermo_thermalExpansionCoeff', tp.tp_id);
        end

        function t = T_sat(tp, p)
            % Get the saturation temperature for a given pressure.
            %
            % parameter p:
            %    Pressure. Unit: Pa.
            % return:
            %    Saturation temperature for pressure p. Unit: K.

            checklib;
            t = calllib(ct, 'thermo_satTemperature', tp.tp_id, p);
        end

        function v = vaporFraction(tp)
            % Get the vapor fractions.
            %
            % return:
            %    Vapor fraction.

            checklib;
            v = calllib(ct, 'thermo_vaporFraction', tp.tp_id);
        end

        function enthalpy = get.H(tp)
            % Get the enthalpy.
            %
            % return:
            %    Enthalpy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                enthalpy = calllib(ct, 'thermo_enthalpy_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                enthalpy = calllib(ct, 'thermo_enthalpy_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function enthalpy = enthalpies_RT(tp)
            % Get the non-dimensional enthalpy.
            %
            % return:
            %    Vector of standard-state species enthalpies divided by RT.

            checklib;
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'thermo_getEnthalpies_RT', tp.tp_id, nsp, pt);
            enthalpy = pt.Value;
        end

        function entropy = get.S(tp)
            % Get the entropy.
            %
            % return:
            %    Entropy of the mixture depending on the basis.
            %    Units: J/kmol-K (molar) J/kg-K (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                entropy = calllib(ct, 'thermo_entropy_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                entropy = calllib(ct, 'thermo_entropy_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function intEnergy = get.U(tp)
            % Get the internal energy.
            %
            % return:
            %    Internal energy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                intEnergy = calllib(ct, 'thermo_intEnergy_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                intEnergy = calllib(ct, 'thermo_intEnergy_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function gibbs = get.G(tp)
            % Get the Gibss free energy.
            %
            % return:
            %    Gibbs free energy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                gibbs = calllib(ct, 'thermo_gibbs_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                gibbs = calllib(ct, 'thermo_gibbs_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        %% ThermoGet multi methods

        function output = get.DP(tp)
            output = {tp.D, tp.P};
        end

        function output = get.DPX(tp)
            % Get density, pressure, and mole fractions depending on the basis.
            % return:
            %    Density. Unit: kmol/m^3 (molar) kg/m^3 (mass).
            %    Pressure. Unit: Pa.
            %    Mole fractions of all species.

            output = {tp.D, tp.P, tp.X};
        end

        function output = get.DPY(tp)
            output = {tp.D, tp.P, tp.Y};
        end

        function output = get.HP(tp)
            output = {tp.H, tp.P};
        end

        function output = get.HPX(tp)
            % Get enthalpy, pressure, and mole fractions depending on the basis.
            % return:
            %    Enthalpy. Unit: J/kmol (molar) J/kg (mass).
            %    Pressure. Unit: Pa.
            %    Mole fractions of all species.

            output = {tp.H, tp.P, tp.X};
        end

        function output = get.HPY(tp)
            output = {tp.H, tp.P, tp.Y};
        end

        function output = get.PV(tp)
            output = {tp.P, tp.V};
        end

        function output = get.PVX(tp)
            output = {tp.P, tp.V, tp.X};
        end

        function output = get.PVY(tp)
            output = {tp.P, tp.V, tp.Y};
        end

        function output = get.SH(tp)
            output = {tp.S, tp.H};
        end

        function output = get.SHX(tp)
            output = {tp.S, tp.H, tp.X};
        end

        function output = get.SHY(tp)
            output = {tp.S, tp.H, tp.Y};
        end

        function output = get.SP(tp)
            output = {tp.S, tp.P};
        end

        function output = get.SPX(tp)
            output = {tp.S, tp.P, tp.X};
        end

        function output = get.SPY(tp)
            output = {tp.S, tp.P, tp.Y};
        end

        function output = get.ST(tp)
            output = {tp.S, tp.T};
        end

        function output = get.STX(tp)
            output = {tp.S, tp.T, tp.X};
        end

        function output = get.STY(tp)
            output = {tp.S, tp.T, tp.Y};
        end

        function output = get.SV(tp)
            output = {tp.S, tp.V};
        end

        function output = get.SVX(tp)
            output = {tp.S, tp.V, tp.X};
        end

        function output = get.SVY(tp)
            output = {tp.S, tp.V, tp.Y};
        end

        function output = get.TD(tp)
            output = {tp.T, tp.D};
        end

        function output = get.TDX(tp)
            output = {tp.T, tp.D, tp.X};
        end

        function output = get.TDY(tp)
            output = {tp.T, tp.D, tp.Y};
        end

        function output = get.TH(tp)
            output = {tp.T, tp.H};
        end

        function output = get.THX(tp)
            output = {tp.T, tp.H, tp.X};
        end

        function output = get.THY(tp)
            output = {tp.T, tp.H, tp.Y};
        end

        function output = get.TP(tp)
            output = {tp.T, tp.P};
        end

        function output = get.TPX(tp)
            output ={tp.T, tp.P, tp.X};
        end

        function output = get.TPY(tp)
            output = {tp.T, tp.P, tp.Y};
        end

        function output = get.TV(tp)
            output = {tp.T, tp.V};
        end

        function output = get.TVX(tp)
            output ={tp.T, tp.V, tp.X};
        end

        function output = get.TVY(tp)
            output = {tp.T, tp.V, tp.Y};
        end

        function output = get.UV(tp)
            output = {tp.U, tp.V};
        end

        function output = get.UVX(tp)
            output = {tp.U, tp.V, tp.X};
        end

        function output = get.UVY(tp)
            output = {tp.U, tp.V, tp.Y};
        end

        function output = get.UP(tp)
            output = {tp.U, tp.P};
        end

        function output = get.UPX(tp)
            output = {tp.U, tp.P, tp.X};
        end

        function output = get.UPY(tp)
            output = {tp.U, tp.P, tp.Y};
        end

        function output = get.VH(tp)
            output = {tp.V, tp.H};
        end

        function output = get.VHX(tp)
            output = {tp.V, tp.H, tp.X};
        end

        function output = get.VHY(tp)
            output = {tp.V, tp.H, tp.Y};
        end

        %% PhaseSet single methods

        function tp = setElectricPotential(tp, phi)
            % Set the electric potential in V.
            calllib(ct, 'thermo_setElectricPotential', tp.tp_id, phi);
        end

        function tp = setState_Psat(tp, p, q)
            % Set saturated vapor
            checklib;
            calllib(ct, 'thermo_setState_Psat', tp.tp_id, p, q);
        end

        function tp = setState_Tsat(tp, t, q)
            % Set saturated liquid
            checklib;
            calllib(ct, 'thermo_setState_Tsat', tp.tp_id, t, 1 - q);
        end

        function set.T(tp, temperature)
            checklib;
            if temperature <= 0
                error('The temperature must be positive');
            end
            calllib(ct, 'thermo_setTemperature', tp.tp_id, temperature);
        end

        function set.P(tp, pressure)
            checklib;
            if pressure <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_setPressure', tp.tp_id, pressure);
        end

        function set.D(tp, density)
            checklib;
            if density <= 0
                error('The density must be positive');
            end
            calllib(ct, 'thermo_setDensity', tp.tp_id, density);
        end

        function set.X(tp, xx)
            checklib;
            lim = 1e-9;
            if isa(xx, 'double')
                nsp = tp.nSpecies;
                if sum(xx) - 1 <= 1e-9
                    norm = 0;
                else norm = 1;
                end
                calllib(ct, 'thermo_setMoleFractions', tp.tp_id, ...
                        nsp, xx, norm);
            elseif isa(xx, 'char')
                calllib(ct, 'thermo_setMoleFractionsByName', tp.tp_id, xx);
            end
        end

        function set.Y(tp, yy)
            checklib;
            if isa(yy, 'double')
                nsp = tp.nSpecies;
                if sum(yy) -1 <= 1e-9
                    norm = 0;
                else norm = 1;
                end
                calllib(ct, 'thermo_setMassFractions', tp.tp_id, ...
                        nsp, yy, norm);
            elseif isa(yy, 'char')
                calllib(ct, 'thermo_setMassFracitonsByName', tp.tp_id, yy);
            end
        end

        %% PhaseSet multi methods

        function set.DP(tp, input)
            checklib;
            d = input{1};
            p = input{2};
            if d <= 0
                error('The density must be positive');
            end
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_RP', tp.tp_id, [d, p]);
        end

        function set.DPX(tp, input)
            tp.X = input{3};
            tp.DP = input(1:2);
        end

        function set.DPY(tp, input)
            tp.Y = input{3};
            tp.DP = input(1:2);
        end

        function set.HP(tp, input)
            checklib;
            h = input{1};
            p = input{2};
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_HP', tp.tp_id, [h, p]);
        end

        function set.HPX(tp, input)
            tp.X = input{3};
            tp.HP = input(1:2);
        end

        function set.HPY(tp, input)
            tp.Y = input{3};
            tp.HP = input(1:2);
        end

        function set.PV(tp, input)
            checklib;
            p = input{1};
            v = input{2};
            if p <= 0
                error('The pressure must be positive');
            end
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_PV', tp.tp_id, [p, v]);
        end

        function set.PVX(tp, input)
            tp.X = input{3};
            tp.PV = input(1:2);
        end

        function set.PVY(tp, input)
            tp.Y = input{3};
            tp.PV = input(1:2);
        end

        function set.SH(tp, input)
            checklib;
            s = input{1};
            h = input{2};
            calllib(ct, 'thermo_set_SH', tp.tp_id, [s, h]);
        end

        function set.SHX(tp, input)
            tp.X = input{3};
            tp.SH = input(1:2);
        end

        function set.SHY(tp, input)
            tp.Y = input{3};
            tp.SH = input(1:2);
        end

        function set.SP(tp, input)
            checklib;
            s = input{1};
            p = input{2};
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_SP', tp.tp_id, [s, p]);
        end

        function set.SPX(tp, input)
            tp.X = input{3};
            tp.SP = input(1:2);
        end

        function set.SPY(tp, input)
            tp.Y = input{3};
            tp.SP = input(1:2);
        end

        function set.ST(tp, input)
            checklib;
            s = input{1};
            t = input{2};
            if t <= 0
                error('The temperature must be positive');
            end
            calllib(ct, 'thermo_set_ST', tp.tp_id, [s, t]);
        end

        function set.STX(tp, input)
            tp.X = input{3};
            tp.ST = input(1:2);
        end

        function set.STY(tp, input)
            tp.Y = input{3};
            tp.ST = input(1:2);
        end

        function set.SV(tp, input)
            checklib;
            s = input{1};
            v = input{2};
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_SV', tp.tp_id, [s, v]);
        end

        function set.SVX(tp, input)
            tp.X = input{3};
            tp.SV = input(1:2);
        end

        function set.SVY(tp, input)
            tp.Y = input{3};
            tp.SP = input(1:2);
        end

        function set.TD(tp, input)
            t = input{1};
            d = input{2};
            if t <= 0
                error('The temperature must be positive');
            end
            if d <= 0
                error('The density must be positive');
            end
            tp.T = t;
            tp.D = d;
        end

        function set.TDX(tp, input)
            tp.X = input{3};
            tp.TD = input(1:2);
        end

        function set.TDY(tp, input)
            tp.Y = input{3};
            tp.TD = input(1:2);
        end

        function set.TH(tp, input)
            checklib;
            t = input{1};
            if t <= 0
                error('The temperature must be positive');
            end
            h = input{2};
            calllib(ct, 'thermo_set_TH', tp.tp_id, [t, h]);
        end

        function set.THX(tp, input)
            tp.X = input{3};
            tp.TH = input(1:2);
        end

        function set.THY(tp, input)
            tp.Y = input{3};
            tp.TH = input(1:2);
        end

        function set.TP(tp, input)
            t = input{1};
            p = input{2};
            if t <= 0
                error('The temperature must be positive');
            end
            if p <= 0
                error('The pressure must be positive');
            end
            tp.T = t;
            tp.P = p;
        end

        function set.TPX(tp, input)
            tp.X = input{3};
            tp.TP = input(1:2);
        end

        function set.TPY(tp, input)
            tp.Y = input{3};
            tp.TP = input(1:2);
        end

        function set.TV(tp, input)
            checklib;
            t = input{1};
            v = input{2};
            if t <= 0
                error('The temperature must be positive');
            end
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_TV', tp.tp_id, [t, v]);
        end

        function set.TVX(tp, input)
            tp.X = input{3};
            tp.TV = input(1:2);
        end

        function set.TVY(tp, input)
            tp.Y = input{3};
            tp.TV = input(1:2);
        end

        function set.UP(tp, input)
            checklib;
            u = input{1};
            p = input{2};
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_UP', tp.tp_id, [u, p]);
        end

        function set.UPX(tp, input)
            tp.X = input{3};
            tp.UP = input(1:2);
        end

        function set.UPY(tp, input)
            tp.Y = input{3};
            tp.UP = input(1:2);
        end

        function set.UV(tp, input)
            checklib;
            u = input{1};
            v = input{2};
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_UV', tp.tp_id, [u, v]);
        end

        function set.UVX(tp, input)
            tp.X = input{3};
            tp.UV = input(1:2);
        end

        function set.UVY(tp, input)
            tp.Y = input{3};
            tp.UV = input(1:2);
        end

        function set.VH(tp, input)
            checklib;
            v = input{1};
            h = input{2};
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_VH', tp.tp_id, [v, h]);
        end

        function set.VHX(tp, input)
            tp.X = input{3};
            tp.VH = input(1:2);
        end

        function set.VHY(tp, input)
            tp.Y = input{3};
            tp.VH = input(1:2);
        end

        function tp = equilibrate(tp, xy, solver, rtol, maxsteps, ...
                                  maxiter, loglevel)
            checklib;
            % use the ChemEquil solver by default
            if nargin < 3
                solver = -1;
            end
            if nargin < 4
                rtol = 1.0e-9;
            end
            if nargin < 5
                maxsteps = 1000;
            end
            if nargin < 6
                maxiter = 100;
            end
            if nargin < 7
                loglevel = 0;
            end
            calllib(ct, 'thermo_equilibrate', tp.tp_id, xy, solver, ...
                    rtol, maxsteps, maxiter, loglevel);
        end

    end
end

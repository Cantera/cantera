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
%         V % specific volume not added yet
        basis
    end

    properties (Constant = true)
        lib = 'cantera_shared'
    end

    properties (Dependent)
        DP
        DPX
        DPY
        HP
        HPX
        HPY
        SP
        SPX
        SPY
%         SVX
%         SVY
        TD
        TDX
        TDY
        TP
        TPX
        TPY
%         UV
%         UVX
%         UVY
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
            tp.tp_id = calllib(tp.lib, 'thermo_newFromFile', src, id);
            tp.basis = 'molar';
        end

        %% Utility methods

        function display(tp, threshold)
            % Display thermo properties

            if nargin < 2 || ~isnumeric(threshold)
                threshold = 1e-14;
            end
            calllib(tp.lib, 'thermo_print', tp.tp_id, 1, threshold);
        end

        function clear(tp)
            % Delete the kernel object.

            checklib;
            calllib(tp.lib, 'thermo_del', tp.tp_id);
        end

        function tp = set.basis(tp, b)
            % Set basis of thermodynamic properties
            % Default is molar

            if strcmp(b, 'mole') || strcmp(b, 'molar') ...
               || strcmp(b, 'Mole') || strcmp(b, 'Molar')
                tp.basis = 'molar'
            elseif strcmp(b, 'mass') || strcmp(b, 'Mass')
                tp.basis = 'mass'
            else error("Basis must be mass or molar")
            end
        end

        %% PhaseGet single methods

        function temperature = get.T(tp)
            % Get the temperature.
            %
            % :return:
            %    Double temperature. Unit: K

            checklib;
            temperature = calllib(tp.lib, 'thermo_temperature', tp.tp_id);
        end

        function pressure = get.P(tp)
            % Get the pressure.
            %
            % :return:
            %    Double pressure. Unit: Pa

            checklib;
            pressure = calllib(tp.lib, 'thermo_pressure', tp.tp_id);
        end

        function density = get.D(tp)
            % Get the density depending on the basis.
            %
            % :return:
            %    Density depending on the basis. Units:
            %    kmol/m^3 (molar) kg/m^3 (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                density = calllib(tp.lib, 'thermo_density', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                density = calllib(tp.lib, 'thermo_molarDensity', tp.tp_id);
            else error("basis not specified");
            end
        end

        function mmw = meanMolecularWeight(tp)
            % Get the mean molecular weight.
            %
            % :return:
            %    Double mean molecular weight. Unit: kg/kmol

            checklib;
            mmw = calllib(tp.lib, 'thermo_meanMolecularWeight', tp.tp_id);
        end

        function nel = nElements(tp)
            % Get the number of elements.
            %
            % :return:
            %    Integer number of elements in the phase.

            checklib;
            nel = calllib(tp.lib, 'thermo_nElements', tp.tp_id);
        end

        function nsp = nSpecies(tp)
            % Get the number of species.
            %
            % :return:
            %    Integer number of species in the phase.

            checklib;
            nsp = calllib(tp.lib, 'thermo_nSpecies', tp.tp_id);
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
            % :param name:
            %    String or cell array of species whose index is requested.
            % :return:
            %    Integer number of species in the phase.

            checklib;
            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);
                for i = 1:m
                    for j = 1:n
                        k(i, j) = calllib(tp.lib, 'thermo_speciesIndex', ...
                                          tp.tp_id, name{i, j}) + 1;
                        if k(i, j) > 1e3
                            warning(['Species ', name{i, j}, ...
                                   ' does not exist in the phase']);
                            k(i, j) = -1;
                        end
                    end
                end
            elseif ischar(name)
                k = calllib(tp.lib, 'thermo_speciesIndex', ...
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

        function k = elementIndex(tp, name)
            % Get the index of an element given the name.
            % The index is an integer assigned to each element in sequence
            % as it is read in from the input file.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.
            %
            % :param name:
            %    String or cell array of elements whose index is requested
            % :return:
            %    Integer number of elements in the phase.

            checklib;
            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);
                for i = 1:m
                    for j = 1:n
                        k(i, j) = calllib(tp.lib, 'thermo_elementIndex', ...
                                          tp.tp_id, name{i, j}) + 1;
                        if k(i, j) > 1e3
                            warning(['Element ', name{i, j}, ...
                                   ' does not exist in the phase']);
                            k(i, j) = -1;
                        end
                    end
                end
            elseif ischar(name)
                k = calllib(tp.lib, 'thermo_elementIndex', ...
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

        function n = nAtoms(tp, species, element)
            % Get the number of atoms of an element in a species.
            %
            % :param k:
            %    String species name or integer species number.
            % :param m:
            %    String element name or integer element number.
            % :return:
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
                n = calllib(tp.lib, 'thermo_nAtoms', tp.tp_id, k-1, m-1);
            else
                error('Two input arguments required.')
            end
        end

        function moleFractions = get.X(tp)
            % Get the mole fractions of all species.
            %
            % :return:
            %    Vector of species mole fractions.

            checklib;
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(tp.lib, 'thermo_getMoleFractions', ...
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
            % :param species:
            %    String or cell array of species whose mole fraction is
            %    requested.
            % :return:
            %    Scalar or vector of species mole fractions.

            xarray = tp.X;
            if isa(species, 'char')
                k = tp.speciesIndex(tp, species);
                if  k > 0
                    x = xarray(k);
                else error("species not found.");
                end
            elseif isa(species, 'cell')
                n = length(species);
                x = zeros(1, n);
                for j = 1:n
                    k = tp.speciesIndex(tp, species{j});
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
            % :return:
            %   Vector of species mass fractions.

            checklib;
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(tp.lib, 'thermo_getMassFractions', ...
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
            % :param species:
            %    String or cell array of species whose mass fraction is
            %    requested.
            % :return:
            %    Scalar or vector of species mass fractions.

            yarray = tp.Y;
            if isa(species, 'char')
                k = tp.speciesIndex(tp, species);
                if  k > 0
                    y = yarray(k);
                else error("Error: species not found.")
                end
            elseif isa(species, 'cell')
                n = length(species);
                y = zeros(1, n);
                for j = 1:n
                    k = tp.speciesIndex(tp, species{j});
                    if k > 0
                        y(j) = yarray(k);
                    else error("Error: species not found.")
                    end
                end
            end
        end

        function mw = MolecularWeights(tp)
            % Get the array of molecular weights of all species.
            %
            % :return:
            %   Vector of species molecular weights. Unit: kg/kmol

            checklib;
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(tp.lib, 'thermo_getMolecularWeights', ...
                    tp.tp_id, nsp, pt);
            mw = pt.Value;
        end

        function e = charges(tp)
            % Get the array of species charges.
            %
            % :return:
            %    Vector of species charges. Unit: elem. charge

            checklib;
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(tp.lib, 'thermo_getCharges', ...
                    tp.tp_id, nsp, pt);
            e = pt.Value;
        end

        function amu = atomicMasses(tp)
            % Get the atomic masses of the elements.
            %
            % :return:
            %    Vector of element atomic masses. Unit: kg/kmol

            checklib;
            nel = tp.nElements;
            aa = zeros(1, nel);
            pt = libpointer('doublePtr', aa);
            calllib(tp.lib, 'thermo_getAtomicWeights', ...
                    tp.tp_id, nel, pt);
            amu = pt.Value;
        end

        function nm = speciesName(tp, k) %need to solve string ptr issue
            % Get the name of a species given the index.
            %
            % :param k:
            %    Scalar or array of integer species index.
            % :return:
            %    Cell array of strings species name.
            [m, n] = size(k);
            nm = cell(m, n);
            for i = 1:m
                for j = 1:n
                    ksp = k(i, j) - 1;
                    buflen = calllib(tp.lib, 'thermo_getSpeciesName', ...
                                     tp.tp_id, ksp, 0, '');
                    if buflen > 0
                        aa = char(zeros(1, buflen));
                        pt = libpointer('voidPtr', int8(aa));
                        pt2 = libpointer('cstring', aa);
                        out_buf = calllib(tp.lib, 'thermo_getSpeciesName', ...
                                          tp.tp_id, ksp, buflen, pt);
                        out2_buf = calllib(tp.lib, 'thermo_getSpeciesName', ...
                                           tp.tp_id, ksp, buflen, pt2);
                        out3_buf = calllib(tp.lib, 'thermo_getSpeciesName', ...
                                           tp.tp_id, ksp, buflen, aa);
                        nm = pt2.Value;
                    end
                end
            end
        end

        %% ThermoGet single methods

        function c = cv(tp)
            % Get the specific heat at constant volume.
            %
            % :return:
            %    Specific heat of the mixture at constant volume depending
            %    on the basis. Units: J/kmol-K (molar) J/kg-K (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                c = calllib(tp.lib, 'thermo_cv_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                c = calllib(tp.lib, 'thermo_cv_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function c = cp(tp)
            % Get the specific heat at constant pressure.
            %
            % :return:
            %    Specific heat of the mixture at constant pressure depending
            %    on the basis. Units: J/kmol-K (molar) J/kg-K (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                c = calllib(tp.lib, 'thermo_cp_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                c = calllib(tp.lib, 'thermo_cp_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function enthalpy = get.H(tp)
            % Get the enthalpy.
            %
            % :return:
            %    Enthalpy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                enthalpy = calllib(tp.lib, 'thermo_enthalpy_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                enthalpy = calllib(tp.lib, 'thermo_enthalpy_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function entropy = get.S(tp)
            % Get the entropy.
            %
            % :return:
            %    Entropy of the mixture depending on the basis.
            %    Units: J/kmol-K (molar) J/kg-K (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                entropy = calllib(tp.lib, 'thermo_entropy_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                entropy = calllib(tp.lib, 'thermo_entropy_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function intEnergy = get.U(tp)
            % Get the internal energy.
            %
            % :return:
            %    Internal energy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                intEnergy = calllib(tp.lib, 'thermo_intEnergy_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                intEnergy = calllib(tp.lib, 'thermo_intEnergy_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function gibbs = G(tp)
            % Get the Gibss free energy.
            %
            % :return:
            %    Gibbs free energy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            checklib;
            if strcmp(tp.basis, 'molar')
                gibbs = calllib(tp.lib, 'thermo_gibbs_mole', tp.tp_id);
            elseif strcmp(tp.basis, 'mass')
                gibbs = calllib(tp.lib, 'thermo_gibbs_mass', tp.tp_id);
            else error("basis not specified");
            end
        end

        function t = minTemp(tp)
            % Get the minimum temperature of the parameter fits.
            %
            % :return:
            %    Vector of minimum temperatures of all species.

            checklib;
            t = calllib(tp.lib, 'thermo_minTemp', tp.tp_id, -1);
        end

        function t = maxTemp(tp)
            % Get the maximum temperature of the parameter fits.
            %
            % :return:
            %    Vector of maximum temperatures of all species.

            checklib;
            t = calllib(tp.lib, 'thermo_maxTemp', tp.tp_id, -1);
        end

        function p = refPressure(tp)
            % Get the reference pressure.
            %
            % :return:
            %    Reference pressure. Unit: Pa.

            checklib;
            p = calllib(tp.lib, 'thermo_refPressure', tp.tp_id, -1);
        end

        function t = critTemperature(tp)
            % Get the critical temperature.
            %
            % :return:
            %    Critical temperature. Unit: K.
            checklib;
            t = calllib(tp.lib, 'thermo_critTemperature', tp.tp_id);
         end

        function p = critPressure(tp)
            % Get the critical pressure.
            %
            % :return:
            %    Critical temperature. Unit: Pa.

            checklib;
            p = calllib(tp.lib, 'thermo_critPressure', tp.tp_id);
        end

        function d = critDensity(tp)
            % Get the critical density.
            %
            % :return:
            %    Critical density. Unit: K.

            checklib;
            d = calllib(tp.lib, 'thermo_critDensity', tp.tp_id);
        end

        function v = vaporFraction(tp)
            % Get the vapor fractions.
            %
            % :return:
            %    Vapor fraction.

            checklib;
            v = calllib(tp.lib, 'thermo_vaporFraction', tp.tp_id);
        end

        function t = T_sat(tp, p)
            % Get the saturation temperature for a given pressure.
            %
            % :param p:
            %    Pressure. Unit: Pa.
            % :return:
            %    Saturation temperature for pressure p. Unit: K.

            checklib;
            t = calllib(tp.lib, 'thermo_satTemperature', tp.tp_id, p);
        end

        function p = P_sat(tp, t)
            % Get the saturation pressure for a given temperature.
            %
            % :param t:
            %    Temperature. Unit: K.
            % :return:
            %    Saturation pressure for temperature t. Unit: Pa.

            checklib;
            p = calllib(tp.lib, 'thermo_satPressure', tp.tp_id, t);
        end

        function v = electricPotential(tp)
            % Get the electric potential
            %
            % :return:
            %    Electric potential of the phase. Unit: V.

            checklib;
            v = calllib(tp.lib, 'thermo_electricPotential', tp.tp_id);
        end

        function b = isothermalCompressibility(tp)
            % Get the isothermal compressibility
            %
            % :return:
            %    Isothermal compressibility. Unit: 1/Pa.

            checklib;
            b = calllib(tp.lib, 'thermo_isothermalCompressibility', tp.tp_id);
        end

        function a = thermalExpansionCoeff(tp)
            % Get the thermal expansion coefficient.
            %
            % :return:
            %    Thermal expansion coefficient. Unit: 1/K.

            checklib;
            a = calllib(tp.lib, 'thermo_thermalExpansionCoeff', tp.tp_id);
        end

        function mu = chemical_potentials(tp)
            % Get the chemical potentials of the species.
            %
            % :return:
            %    Vector of species chemical potentials. Unit: J/kmol.

            checklib;
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(tp.lib, 'thermo_chemPotentials', ...
                    tp.tp_id, nsp, pt);
            moleFractions = pt.Value;
        end

        %% ThermoGet multi methods

        function output = get.DP(tp)
            % Get density and pressure depending on the basis.
            % :return:
            %    Density. Unit: kmol/m^3 (molar) kg/m^3 (mass)
            %    Pressure. Unit: Pa

            output = [tp.D, tp.P];
        end

        function output = get.DPX(tp)
            % Get density, pressure, and mole fractions depending on the basis.
            % :return:
            %    Density. Unit: kmol/m^3 (molar) kg/m^3 (mass).
            %    Pressure. Unit: Pa.
            %    Mole fractions of all species.

            output = [tp.D, tp.P, tp.X];
        end

        function output = get.DPY(tp)
            % Get density, pressure, and mass fractions depending on the basis.
            % :return:
            %    Density. Unit: kmol/m^3 (molar) kg/m^3 (mass).
            %    Pressure. Unit: Pa.
            %    Mass fractions of all species.

            output = [tp.D, tp.P, tp.Y];
        end

        function output = get.HP(tp)
            % Get enthalpy and pressure depending on the basis.
            % :return:
            %    Enthalpy. Unit: J/kmol (molar) J/kg (mass).
            %    Pressure. Unit: Pa.

            output = [tp.H, tp.P];
        end

        function output = get.HPX(tp)
            % Get enthalpy, pressure, and mole fractions depending on the basis.
            % :return:
            %    Enthalpy. Unit: J/kmol (molar) J/kg (mass).
            %    Pressure. Unit: Pa.
            %    Mole fractions of all species.

            output = [tp.H, tp.P, tp.X];
        end

        function output = get.HPY(tp)
            % Get enthalpy, pressure, and mole fractions depending on the basis.
            % :return:
            %    Enthalpy. Unit: J/kmol (molar) J/kg (mass).
            %    Pressure. Unit: Pa.
            %    Mass fractions of all species.

            output = [tp.H, tp.P, tp.Y];
        end

        function output = get.SP(tp)
            % Get entropy, pressure, and mole fractions depending on the basis.
            % :return:
            %    Entropy. Unit: J/kmol-K (molar) J/kg-K (mass).
            %    Pressure. Unit: Pa.

            output = [tp.S, tp.P];
        end

        function output = get.SPX(tp)
            % Get entropy, pressure, and mole fractions depending on the basis.
            % :return:
            %    Entropy. Unit: J/kmol-K (molar) J/kg-K (mass).
            %    Pressure. Unit: Pa.
            %    Mole fractions of all species.

            output = [tp.S, tp.P, tp.X];
        end

        function output = get.SPY(tp)
            % Get entropy, pressure, and mass fractions depending on the basis.
            % :return:
            %    Entropy. Unit: J/kmol-K (molar) J/kg-K (mass).
            %    Pressure. Unit: Pa.
            %    Mass fractions of all species.

            output = [tp.H, tp.P, tp.Y];
        end

        function output = get.TD(tp)
            % Get temperature and density depending on the basis.
            % :return:
            %    Temperature. Unit: K.
            %    Density. Unit: kmol/m^3 (molar) kg/m^3 (mass).

            output = [tp.T, tp.D];
        end

        function output = get.TDX(tp)
            % Get entropy, pressure, and mole fractions depending on the basis.
            % :return:
            %    Temperature. Unit: K.
            %    Density. Unit: kmol/m^3 (molar) kg/m^3 (mass).
            %    Mole fractions of all species.

            output = [tp.T, tp.D, tp.X];
        end

        function output = get.TDY(tp)
            % Get entropy, pressure, and mass fractions depending on the basis.
            % :return:
            %    Temperature. Unit: K.
            %    Density. Unit: kmol/m^3 (molar) kg/m^3 (mass).
            %    Mass fractions of all species.

            output = [tp.T, tp.D, tp.Y];
        end

        function output = get.TP(tp)
            % Get temperature and density depending on the basis.
            % :return:
            %    Temperature. Unit: K.
            %    Pressure. Unit: Pa.
            output = [tp.T, tp.P];
        end

        function output = get.TPX(tp)
            % Get entropy, pressure, and mole fractions depending on the basis.
            % :return:
            %    Temperature. Unit: K.
            %    Pressure. Unit: Pa.
            %    Mole fractions of all species.

            output = [tp.T, tp.P, tp.X];
        end

        function output = get.TPY(tp)
            % Get entropy, pressure, and mass fractions depending on the basis.
            % :return:
            %    Temperature. Unit: K.
            %    Pressure. Unit: Pa.
            %    Mass fractions of all species.

            output = [tp.T, tp.P, tp.Y];
        end

        %% PhaseSet Single Methods

        function tp = set.T(tp, temperature)
            checklib;
            if temperature <= 0
                error('The temperature must be positive');
            end
            calllib(tp.lib, 'thermo_setTemperature', tp.tp_id, temperature);
        end

        function tp = set.P(tp, pressure)
            checklib;
            if pressure <= 0
                error('The pressure must be positive');
            end
            calllib(tp.lib, 'thermo_setPressure', tp.tp_id, pressure);
        end

        function tp = set.D(tp, density)
            checklib;
            if density <= 0
                error('The density must be positive');
            end
            calllib(tp.lib, 'thermo_setDensity', tp.tp_id, density);
        end

%         function tp = set.X(tp,

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
            calllib(tp.lib, 'thermo_equilibrate', tp.tp_id, xy, solver, ...
                    rtol, maxsteps, maxiter, loglevel);
        end

    end
end

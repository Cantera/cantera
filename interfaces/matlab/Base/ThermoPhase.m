classdef (Abstract) ThermoPhase < handle
    % ThermoPhase Class ::
    %
    %     >> t = ThermoPhase(id)
    %
    % Retrieve instance of class :mat:class:`ThermoPhase` associated with a
    % :mat:class:`Solution` object. The constructor is called whenever a new
    % :mat:class:`Solution` is instantiated and should not be used directly.
    %
    % :param id:
    %     Integer ID of the solution holding the :mat:class:`ThermoPhase` object.

    properties (SetAccess = public)

        X % Mole fractions.

        Y % Mass fractions.

        % Basis ::
        %
        %     >> tp.basis = b
        %
        % Determines whether intensive thermodynamic properties are
        % treated on a mass (per kg) or molar (per kmol) basis. This
        % affects the values returned by the properties H, U, S, G, V,
        % Density, Cv, and Cp, as well as the values used with the
        % state-setting properties such as HPX and UV.
        %
        % :param b:
        %     String. Can be ``mole`` / ``molar`` / ``Molar`` / ``Mole`` or
        %     ``mass`` / ``Mass``.
        basis

        name % Name of the phase.

        electricPotential % Electric potential [V].

    end

    properties (SetAccess = protected)

        tpID % ID of the ThermoPhase object.

    end

    properties (SetAccess = immutable)

        atomicWeights % Atomic weights of the elements [kg/kmol].

        charges % Species charges in units of the elementary charge.

        % Mean molecular weight [kg/kmol] of the mixture
        %
        % The mean molecular weight is the mole-fraction-weighted sum of the
        % atomic weights of the individual species in the phase.
        meanMolecularWeight

        massDensity % Mass basis density [kg/m³].

        molarDensity % Molar basis density [kmol/m³].

        molecularWeights % Molecular weights [kg/kmol] of the species.

        nElements % Number of elements in the phase.

        nSpecies % Number of species in the phase.

        T % Temperature [K].

        P % Pressure [Pa].

        D % Density depending on the basis [kmol/m³ or kg/m³].

        H % Enthalpy depending on the basis [J/kmol or J/kg].

        S % Entropy depending on the basis [J/kmol/K or J/kg/K].

        U % Internal energy depending on the basis [J/kmol or J/kg].

        G % Gibbs free energy depending on the basis [J/kmol or J/kg].

        Q % Vapor fraction of the phase.

        % Basis-dependent specific heat at constant volume and composition
        % [J/kmol/K or J/kg/K].
        cv

        % Basis-dependent specific heat at constant pressure and composition
        % [J/kmol/K or J/kg/K].
        cp

        % Concentrations of the species [kmol/m³ for bulk phases; kmol/m² for surface
        % phases].
        concentrations

        % Chemical potentials of the species [J/kmol].
        chemicalPotentials

        % Electrochemical potentials of the species [J/kmol].
        electrochemicalPotentials

        % Partial molar enthalpies for the species in the mixture [J/kmol].
        partialMolarEnthalpies

        % Partial molar entropies for the species in the mixture [J/kmol/K].
        partialMolarEntropies

        % Partial molar internal energies for the species in the mixture [J/kmol].
        partialMolarIntEnergies

        % Partial molar heat capacities for the species in the mixture [J/kmol/K].
        partialMolarCp

        % Partial molar volumes for the species in the mixture [m³/kmol].
        partialMolarVolumes

        critDensity % Critical density [kg/m³].

        critPressure % Critical pressure [Pa].

        critTemperature % Critical temperature [K].

        eosType % Type of equation of state.

        isIdealGas % A flag indicating whether the phase is an ideal gas.

        isothermalCompressibility % Isothermal compressibility [1/Pa].

        % Maximum temperature for which thermodynamic parameter fits are valid
        % for all species.
        %
        % The parameterizations used to represent the temperature-dependent
        % species thermodynamic properties are generally only valid in some
        % finite temperature range, which may be different for each species
        % in the phase.
        %
        % See also: :mat:meth:`minTemp`
        maxTemp

        % Minimum temperature for which thermodynamic parameter fits are valid
        % for all species.
        %
        % The parameterizations used to represent the temperature-dependent
        % species thermodynamic properties are generally only valid in some
        % finite temperature range, which may be different for each species
        % in the phase.
        %
        % See also: :mat:class:`maxTemp`
        minTemp

        refPressure % Reference pressure [Pa] for standard-state.

        satPressure % Saturation pressure [Pa] at the current temperature

        satTemperature % Saturation temperature [K] at the current pressure

        % Speed of sound [m/s] ::
        %
        %     >> c = tp.soundspeed
        %
        % If the phase is an ideal gas, the speed of sound is calculated by:
        %
        % .. math:: c = \sqrt{\gamma * R * T}
        %
        % where :math:`\gamma` is the ratio of specific heats, :math:`R` is
        % the specific gas constant, and :math:`T` is the temperature. If the
        % phase is not an ideal gas, the speed of sound is calculated by
        %
        % .. math:: c = \sqrt{\left(\frac{\partial p}{\partial \rho}\right)_s}
        %
        % where :math:`p` is the pressure and :math:`\rho` is the density,
        % and the subscript :math:`s` indicates constant entropy. This is
        % approximated by slightly increasing the density at constant entropy
        % and computing the change in pressure.
        %
        % .. math:: c = \sqrt{\frac{p_1 - p_0}{\rho_1-\rho_0}}
        soundSpeed

        speciesNames % Cell array of species names

        thermalExpansionCoeff % Thermal expansion coefficient [1/K].

    end

    properties (Dependent = true)

        V % Basis-dependent specific volume [m³/kmol or m³/kg]

        % Get/Set density [kg/m³ or kmol/m³] and pressure [Pa].
        DP

        % Get/Set density [kg/m³ or kmol/m³], pressure [Pa], and mole fractions.
        DPX

        % Get/Set density [kg/m³ or kmol/m³], pressure [Pa], and mass fractions.
        DPY

        % Get density [kg/m³ or kmol/m³], pressure [Pa], and vapor fraction.
        DPQ

        % Get/Set enthalpy [J/kg or J/kmol] and pressure [Pa].
        HP

        % Get/Set enthalpy [J/kg or J/kmol], pressure [Pa], and mole fractions.
        HPX

        % Get/Set enthalpy [J/kg or J/kmol], pressure [Pa], and mass fractions.
        HPY

        % Get enthalpy [J/kg or J/kmol], pressure [Pa], and vapor fraction.
        HPQ

        % Get/Set pressure [Pa] and specific volume [m³/kg or m³/kmol].
        PV

        % Get/Set pressure [Pa], specific volume [m³/kg or m³/kmol],
        % and mole fractions.
        PVX

        % Get/Set pressure [Pa], specific volume [m³/kg or m³/kmol],
        % and mass fractions.
        PVY

        % Get/Set pressure [Pa] and vapor fraction of a two-phase state.
        PQ

        % Get/Set entropy [J/kg/K or J/kmol/K] and enthalpy [J/kg or J/kmol].
        SH

        % Get/Set entropy [J/kg/K or J/kmol/K], enthalpy [J/kg or J/kmol],
        % and mole fractions.
        SHX

        % Get/Set entropy [J/kg/K or J/kmol/K], enthalpy [J/kg or J/kmol],
        % and mass fractions.
        SHY

        % Get/Set entropy [J/kg/K or J/kmol/K] and pressure [Pa].
        SP

        % Get/Set entropy [J/kg/K or J/kmol/K], pressure [Pa], and mole fractions.
        SPX

        % Get/Set entropy [J/kg/K or J/kmol/K], pressure [Pa], and mass fractions.
        SPY

        % Get entropy [J/kg/K or J/kmol/K], pressure [Pa], and vapor fraction.
        SPQ

        % Get/Set entropy [J/kg/K or J/kmol/K] and temperature [K].
        ST

        % Get/Set entropy [J/kg/K or J/kmol/K], temperature [K], and mole fractions.
        STX

        % Get/Set entropy [J/kg/K or J/kmol/K], temperature [K], and mass fractions.
        STY

        % Get/Set entropy [J/kg/K or J/kmol/K] and specific volume [m³/kg or m³/kmol].
        SV

        % Get/Set entropy [J/kg/K or J/kmol/K], specific volume [m³/kg or m³/kmol],
        % and mole fractions.
        SVX

        % Get/Set entropy [J/kg/K or J/kmol/K], specific volume [m³/kg or m³/kmol],
        % and mass fractions.
        SVY

        % Get/Set entropy [J/kg/K or J/kmol/K], specific volume [m³/kg or m³/kmol],
        % and vapor fraction.
        SVQ

        % Get/Set temperature [K] and density [kg/m³ or kmol/m³].
        TD

        % Get/Set temperature [K], density [kg/m³ or kmol/m³], and mole fractions.
        TDX

        % Get/Set temperature [K], density [kg/m³ or kmol/m³], and mass fractions.
        TDY

        % Get temperature [K], density [kg/m³ or kmol/m³], and vapor fraction.
        TDQ

        % Get/Set temperature [K] and enthalpy [J/kg or J/kmol].
        TH

        % Get/Set temperature [K], enthalpy [J/kg or J/kmol], and mole fractions.
        THX

        % Get/Set temperature [K], enthalpy [J/kg or J/kmol], and mass fractions.
        THY

        % Get/Set temperature [K] and pressure [Pa].
        TP

        % Get/Set temperature [K], pressure [Pa], and mole fractions.
        TPX

        % Get/Set temperature [K], pressure [Pa], and mass fractions.
        TPY

        % Get temperature [K], pressure [Pa], and vapor fraction.
        TPQ

        % Get/Set temperature [K] and vapor fraction of a two-phase state.
        TQ

        % Get/Set temperature [K] and specific volume [m³/kg or m³/kmol].
        TV

        % Get/Set temperature [K], specific volume [m³/kg or m³/kmol],
        % and mole fractions.
        TVX

        % Get/Set temperature [K], specific volume [m³/kg or m³/kmol],
        % and mass fractions.
        TVY

        % Get/Set internal energy [J/kg or J/kmol] and specific volume
        % [m³/kg or m³/kmol].
        UV

        % Get/Set internal energy [J/kg or J/kmol], specific volume
        % [m³/kg or m³/kmol], and mole fractions.
        UVX

        % Get/Set internal energy [J/kg or J/kmol], specific volume
        % [m³/kg or m³/kmol], and mass fractions.
        UVY

        % Get internal energy [J/kg or J/kmol], specific volume
        % [m³/kg or m³/kmol], and vapor fraction.
        UVQ

        % Get/Set internal energy [J/kg or J/kmol] and pressure [Pa].
        UP

        % Get/Set internal energy [J/kg or J/kmol], pressure [Pa], and mole fractions.
        UPX

        % Get/Set internal energy [J/kg or J/kmol], pressure [Pa], and mass fractions.
        UPY

        % Get/Set volume [m³/kg or m³/kmol] and enthalpy [J/kg or J/kmol].
        VH

        % Get/Set specific volume [m³/kg or m³/kmol], enthalpy [J/kg or J/kmol],
        % and mole fractions.
        VHX

        % Get/Set specific volume [m³/kg or m³/kmol], enthalpy [J/kg or J/kmol],
        % and mass fractions.
        VHY

    end

    methods
        %% ThermoPhase Class Constructor

        function tp = ThermoPhase(id)
            % Create a :mat:class:`ThermoPhase` object.
            arguments
                id (1,1) double {mustBeInteger}
            end

            tp.tpID = ctFunc('mSol_thermo', id);
            tp.basis = 'molar';
        end

        %% ThermoPhase Utility Methods

        function display(obj)
            disp(obj.report);
        end

        function tp = equilibrate(obj, xy, solver, rtol, maxsteps, maxiter, loglevel)
            % Set the phase to a state of chemical equilibrium ::
            %
            %     >> tp.equilibrate(xy, solver, rtol, maxsteps, maxiter, loglevel)
            %
            % :param XY:
            %     A two-letter string, which must be one of the set
            %     ``['TP','TV','HP','SP','SV','UV','UP']``,
            %     indicating which pair of properties should be held constant.
            %     Not all of the properties to be held constant are available with
            %     all of the solvers.
            % :param solver:
            %     Specifies the equilibrium solver to use. Choices are
            %     'element_potential' (fast solver using the element potential method),
            %     'gibbs' (slower but more robust Gibbs minimization solver), 'vcs'
            %     (VCS algorithm). For the default 'auto' setting, the fast solver
            %     will be tried first, then if it fails the Gibbs minimization solver
            %     will be tried.
            % :param rtol:
            %     The relative error tolerance.
            % :param maxsteps:
            %     Maximum number of steps in composition to take to find a
            %     converged solution.
            % :param maxiter:
            %     For the Gibbs minimization solver only, this specifies the number
            %     of 'outer' iterations on T or P when some property pair other than
            %     TP is specified.
            % :param loglevel:
            %     Set to a value > 0 to write diagnostic output. Larger values
            %     generate more detailed information.

            arguments
                obj (1,1) ThermoPhase
                xy (1,1) string {mustBeMember(xy, ["TP", "TV", "HP", "SP", ...
                                                   "SV", "UV", "UP"])} = "TP"
                solver (1,1) string {mustBeMember(solver, ["auto", ...
                                     "vcs", "gibbs", "element_potential"])} = "auto"
                rtol (1,1) double {mustBePositive} = 1.0e-9
                maxsteps (1,1) double {mustBeInteger, mustBePositive} = 1000
                maxiter (1,1) double {mustBeInteger, mustBePositive} = 100
                loglevel (1,1) double {mustBeInteger, mustBeNonnegative} = 0
            end

            ctFunc('mThermo_equilibrate', obj.tpID, xy, solver, rtol, ...
                   maxsteps, maxiter, loglevel);
        end

        %% ThermoPhase inquiry methods

        function k = elementIndex(obj, name)
            % Index of an element given its name ::
            %
            %     >> k = tp.elementIndex(name)
            %
            % The index is an integer assigned to each element in sequence as it
            % is read in from the input file.
            %
            % If ``name`` is a single string, the return value will be a integer
            % containing the corresponding index. If it is an cell array of
            % strings, the output will be an array of the same shape
            % containing the indices.
            %
            % NOTE: In keeping with the conventions used by Matlab, this method
            % returns 1 for the first element. In contrast, the corresponding
            % method elementIndex in the Cantera C++ and Python interfaces
            % returns 0 for the first element, 1 for the second one, etc. ::
            %
            %     >> ic = gas.elementIndex('C');
            %     >> ih = gas.elementIndex('H');
            %
            % :param name:
            %     String or cell array of strings of elements to look up
            % :return:
            %     Integer or vector of integers of element indices

            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);

                for i = 1:m

                    for j = 1:n
                        k(i, j) = ctFunc('mThermo_elementIndex', ...
                                        obj.tpID, name{i, j}) + 1;

                        if k(i, j) > 1e3
                            warning(['Element ', name{i, j}, ...
                                    ' does not exist in the phase']);
                            k(i, j) = -1;
                        end

                    end

                end

            elseif ischar(name)
                k = ctFunc('mThermo_elementIndex', obj.tpID, name) + 1;

                if k > 1e3
                    warning(['Element ', name, ' does not exist in the phase']);
                    k = -1;
                end

            else
                error('name must be either a string or cell array of strings')
            end

        end

        function elMassFrac = elementalMassFraction(obj, element)
            % Elemental mass fraction in gas object ::
            %
            %     >> elMassFrac = tp.elementalMassFraction(element)
            %
            % The elemental mass fraction in a gas object is calculated using
            % the following equation:
            %
            %  .. math:: elMassFrac = \sum_{1}^{nSpecies} \frac{nAtoms(k, m)*Mel(m)*Y(k)}{mw(k)}
            %
            % where :math:`nAtoms(k, m)` is the number of atoms of element :math:`m` in
            % species :math:`k`; :math:`Mel(m)` is the atomic weight of
            % element :math:`m`; :math:`Y(k)` is the mass fraction of
            % species :math:`k`; and :math:`mw(k)` is the molecular weight of
            % species :math:`k`.
            %
            % :param element:
            %     String representing the element name.
            % :return:
            %     Elemental mass fraction within a gas object.
            arguments
                obj
                element
            end
            if ~ischar(element)
                error('Wrong type for element name: must be character.');
            end

            n = obj.nSpecies;
            spec = obj.speciesNames;
            eli = obj.elementIndex(element);
            M = obj.atomicWeights;
            Mel = M(eli);
            MW = obj.molecularWeights;
            natoms = zeros(1, n);
            yy = zeros(1, n);
            % Initialize the element mass fraction as zero.
            elMassFrac = 0.0;
            % Perform summation of elemental mass fraction over all species.
            for i = 1:n
                natoms(i) = obj.nAtoms(spec{i}, element);
                yy(i) = obj.massFraction(spec{i});
                elMassFrac = elMassFrac + (natoms(i) * Mel * yy(i)) / MW(i);
            end

        end

        function n = nAtoms(obj, species, element)
            % Number of atoms of an element in a species ::
            %
            %   >> n = tp.nAtoms(k,m)
            %
            % :param species:
            %     Species name or index
            % :param element:
            %     Element name or index
            % :return:
            %     Number of atoms of the specified element in the species.
            arguments
                obj
                species
                element
            end

            if ischar(species) | isstring(species)
                k = obj.speciesIndex(species);
            else
                k = species;
            end

            if ischar(element) | isstring(element)
                m = obj.elementIndex(element);
            else
                m = element;
            end

            n = ctFunc('mThermo_nAtoms', obj.tpID, k - 1, m - 1);

        end

        function k = speciesIndex(obj, name)
            % Index of a species given the name ::
            %
            %   >> k = tp.speciesIndex(name)
            %
            % The index is an integer assigned to each species in sequence as it
            % is read in from the input file. ::
            %
            %    >> ich4 = gas.speciesIndex('CH4');
            %    >> iho2 = gas.speciesIndex('HO2');
            %
            % .. note::
            %
            %    In keeping with the conventions used by Matlab, this method returns 1
            %    for the first species, 2 for the second, etc. In contrast, the
            %    corresponding method in the Cantera C++ and Python interfaces returns 0
            %    for the first species, 1 for the second one, etc.
            %
            % :param name:
            %     If name is a single string, the return value will be a integer
            %     containing the corresponding index. If it is an cell array of
            %     strings, the output will be an array of the same shape
            %     containing the indices.
            % :return:
            %     Scalar or array of integers

            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);

                for i = 1:m

                    for j = 1:n
                        k(i, j) = ctFunc('mThermo_speciesIndex', ...
                                        obj.tpID, name{i, j}) + 1;

                        if k(i, j) > 1e6
                            warning(['Species ', name{i, j}, ...
                                    ' does not exist in the phase']);
                            k(i, j) = -1;
                        end

                    end

                end

            elseif ischar(name)
                k = ctFunc('mThermo_speciesIndex', obj.tpID, name) + 1;

                if k > 1e6
                    warning(['Species ', name, ' does not exist in the phase.']);
                    k = -1;
                end

            else
                error('Name must be either a string or cell array of strings.')
            end

        end

        function nm = speciesName(obj, k)
            % Name of one or multiple species given the index ::
            %
            %     >> k = tp.speciesName(k)
            %
            % :param k:
            %     Scalar of array of integers of
            %     NOTE: In keeping with the conventions used by Matlab, the indices of
            %     species start with 1 for the first, then 2 for the second, etc.
            % :return:
            %     Cell array of strings of species inquired.

            [m, n] = size(k);
            nm = cell(m, n);

            for i = 1:m

                for j = 1:n
                    ksp = k(i, j) - 1;
                    output = ctString('mThermo_speciesName', obj.tpID, ksp);
                    nm{i, j} = output;
                end

            end

        end

        function x = moleFraction(obj, species)
            % Get the mole fraction of one or a list of species ::
            %
            %     >> x = tp.moleFraction(species)
            %
            % :param species:
            %     String or cell array of strings of species whose mole
            %     fraction is desired
            % :return:
            %     Scalar or vector double mole fractions

            xarray = obj.X;

            if isa(species, 'char')
                k = obj.speciesIndex(species);

                if k > 0
                    x = xarray(k);
                else
                    error("species not found.");
                end

            elseif isa(species, 'cell')
                n = length(species);
                x = zeros(1, n);

                for j = 1:n
                    k = obj.speciesIndex(species{j});

                    if k > 0
                        x(j) = xarray(k);
                    else
                        error("species not found.");
                    end

                end

            else
                error('Species name must be either a string or cell array of strings.')
            end

        end

        function y = massFraction(obj, species)
            % Get the mass fraction of one or a list of species ::
            %
            %     >> y = tp.massFraction(species)
            %
            % :param species:
            %     String or cell array of strings of species whose mass
            %     fraction is desired
            % :return:
            %     Scalar or vector double mass fractions

            yy = obj.Y;

            if isa(species, 'char')
                k = obj.speciesIndex(species);

                if k > 0
                    y = yy(k);
                else
                    error("Error: species not found.")
                end

            elseif isa(species, 'cell')
                n = length(species);
                y = zeros(1, n);

                for j = 1:n
                    k = obj.speciesIndex(species{j});

                    if k > 0
                        y(j) = yy(k);
                    else
                        error("Error: species not found.")
                    end

                end

            else
                error('Species name must be either a string or cell array of strings.')
            end

        end

        function str = report(obj, threshold)
            arguments
                obj
                threshold (1,1) double = 1e-14
            end
            str = ctString('mThermo_report', obj.tpID, 1, threshold);
        end

        %% Single-property getter methods

        function amu = get.atomicWeights(obj)
            nel = obj.nElements;
            amu = ctArray('mThermo_atomicWeights', nel, obj.tpID);
        end

        function e = get.charges(obj)
            nsp = obj.nSpecies;
            e = ctArray('mThermo_getCharges', nsp, obj.tpID);
        end

        function c = get.cv(obj)
            if strcmp(obj.basis, 'molar')
                c = ctFunc('mThermo_cv_mole', obj.tpID);
            else
                c = ctFunc('mThermo_cv_mass', obj.tpID);
            end
        end

        function c = get.cp(obj)
            if strcmp(obj.basis, 'molar')
                c = ctFunc('mThermo_cp_mole', obj.tpID);
            else
                c = ctFunc('mThermo_cp_mass', obj.tpID);
            end
        end

        function d = get.critDensity(obj)
            d = ctFunc('mThermo_critDensity', obj.tpID);
        end

        function p = get.critPressure(obj)
            p = ctFunc('mThermo_critPressure', obj.tpID);
        end

        function t = get.critTemperature(obj)
            t = ctFunc('mThermo_critTemperature', obj.tpID);
        end

        function v = get.electricPotential(obj)
            v = ctFunc('mThermo_electricPotential', obj.tpID);
        end

        function e = get.eosType(obj)
            e = ctString('mThermo_type', obj.tpID);
        end

        function v = get.isIdealGas(obj)
            v = strcmp(obj.eosType, 'ideal-gas');
        end

        function b = get.isothermalCompressibility(obj)
            b = ctFunc('mThermo_isothermalCompressibility', obj.tpID);
        end

        function t = get.maxTemp(obj)
            t = ctFunc('mThermo_maxTemp', obj.tpID, -1);
        end

        function t = get.minTemp(obj)
            t = ctFunc('mThermo_minTemp', obj.tpID, -1);
        end

        function mmw = get.meanMolecularWeight(obj)
            mmw = ctFunc('mThermo_meanMolecularWeight', obj.tpID);
        end

        function density = get.massDensity(obj)
            density = ctFunc('mThermo_density', obj.tpID);
        end

        function density = get.molarDensity(obj)
            density = ctFunc('mThermo_molarDensity', obj.tpID);
        end

        function c = get.concentrations(obj)
            nsp = obj.nSpecies;
            c = ctArray('mThermo_getConcentrations', nsp, obj.tpID);
        end

        function mw = get.molecularWeights(obj)
            nsp = obj.nSpecies;
            mw = ctArray('mThermo_getMolecularWeights', nsp, obj.tpID);
        end

        function nel = get.nElements(obj)
            nel = ctFunc('mThermo_nElements', obj.tpID);
        end

        function nsp = get.nSpecies(obj)
            nsp = ctFunc('mThermo_nSpecies', obj.tpID);
        end

        function p = get.refPressure(obj)
            p = ctFunc('mThermo_refPressure', obj.tpID);
        end

        function p = get.satPressure(obj)
            p = ctFunc('mThermo_satPressure', obj.tpID, obj.T);
        end

        function t = get.satTemperature(obj)
            t = ctFunc('mThermo_satTemperature', obj.tpID, obj.P);
        end

        function c = get.soundSpeed(obj)

            if obj.isIdealGas
                obj.basis = 'mass';
                gamma = obj.cp / obj.cv;
                wtm = obj.meanMolecularWeight;
                r = 8314.4621 / wtm;
                c = sqrt(gamma * r * obj.T);
            else
                rho0 = obj.D;
                p0 = obj.P;
                s0 = obj.S;
                rho1 = 1.001 * rho0;
                obj.SV = {s0, 1/rho0};
                p1 = obj.P;
                dpdrho_s = (p1 - p0) / (rho1 - rho0);
                c = sqrt(dpdrho_s);
            end

        end

        function s = get.name(obj)
            s = ctString('mThermo_name', obj.tpID);
        end

        function n = get.speciesNames(obj)
            n = obj.speciesName(1:obj.nSpecies);
        end

        function a = get.thermalExpansionCoeff(obj)
            a = ctFunc('mThermo_thermalExpansionCoeff', obj.tpID);
        end

        function temperature = get.T(obj)
            temperature = ctFunc('mThermo_temperature', obj.tpID);
        end

        function pressure = get.P(obj)
            pressure = ctFunc('mThermo_pressure', obj.tpID);
        end

        function v = get.Q(obj)
            v = ctFunc('mThermo_vaporFraction', obj.tpID);
        end

        function density = get.D(obj)
            if strcmp(obj.basis, 'mass')
                density = ctFunc('mThermo_density', obj.tpID);
            else
                density = obj.molarDensity;
            end
        end

        function volume = get.V(obj)
            volume = 1 / obj.D;
        end

        function moleFractions = get.X(obj)
            nsp = obj.nSpecies;
            moleFractions = ctArray('mThermo_getMoleFractions', nsp, obj.tpID);
        end

        function massFractions = get.Y(obj)
            nsp = obj.nSpecies;
            massFractions = ctArray('mThermo_getMassFractions', nsp, obj.tpID);
        end

        function enthalpy = get.H(obj)
            if strcmp(obj.basis, 'molar')
                enthalpy = ctFunc('mThermo_enthalpy_mole', obj.tpID);
            else
                enthalpy = ctFunc('mThermo_enthalpy_mass', obj.tpID);
            end
        end

        function mu = get.chemicalPotentials(obj)
            nsp = obj.nSpecies;
            mu = ctArray('mThermo_chemPotentials', nsp, obj.tpID);
        end

        function emu = get.electrochemicalPotentials(obj)
            nsp = obj.nSpecies;
            emu = ctArray('mThermo_electrochemPotentials', nsp, obj.tpID);
        end

        function enthalpies = get.partialMolarEnthalpies(obj)
            nsp = obj.nSpecies;
            enthalpies = ctArray('mThermo_getPartialMolarEnthalpies', nsp, obj.tpID);
        end

        function entropies = get.partialMolarEntropies(obj)
            nsp = obj.nSpecies;
            entropies = ctArray('mThermo_getPartialMolarEntropies', nsp, obj.tpID);
        end

        function intEnergies = get.partialMolarIntEnergies(obj)
            nsp = obj.nSpecies;
            intEnergies = ctArray('mThermo_getPartialMolarIntEnergies', nsp, obj.tpID);
        end

        function cps = get.partialMolarCp(obj)
            nsp = obj.nSpecies;
            cps = ctArray('mThermo_getPartialMolarCp', nsp, obj.tpID);
        end

        function volumes = get.partialMolarVolumes(obj)
            nsp = obj.nSpecies;
            volumes = ctFunc('mThermo_getPartialMolarVolumes', nsp, obj.tpID);
        end

        function entropy = get.S(obj)
            if strcmp(obj.basis, 'molar')
                entropy = ctFunc('mThermo_entropy_mole', obj.tpID);
            else
                entropy = ctFunc('mThermo_entropy_mass', obj.tpID);
            end
        end

        function intEnergy = get.U(obj)
            if strcmp(obj.basis, 'molar')
                intEnergy = ctFunc('mThermo_intEnergy_mole', obj.tpID);
            else
                intEnergy = ctFunc('mThermo_intEnergy_mass', obj.tpID);
            end
        end

        function gibbs = get.G(obj)
            if strcmp(obj.basis, 'molar')
                gibbs = ctFunc('mThermo_gibbs_mole', obj.tpID);
            else
                gibbs = ctFunc('mThermo_gibbs_mass', obj.tpID);
            end
        end

        %% Multi-property getter methods

        function output = get.DP(obj)
            output = {obj.D, obj.P};
        end

        function output = get.DPX(obj)
            output = {obj.D, obj.P, obj.X};
        end

        function output = get.DPY(obj)
            output = {obj.D, obj.P, obj.Y};
        end

        function output = get.DPQ(obj)
            output = {obj.D, obj.P, obj.Q};
        end

        function output = get.HP(obj)
            output = {obj.H, obj.P};
        end

        function output = get.HPX(obj)
            output = {obj.H, obj.P, obj.X};
        end

        function output = get.HPY(obj)
            output = {obj.H, obj.P, obj.Y};
        end

        function output = get.HPQ(obj)
            output = {obj.H, obj.P, obj.Q};
        end

        function output = get.PV(obj)
            output = {obj.P, obj.V};
        end

        function output = get.PVX(obj)
            output = {obj.P, obj.V, obj.X};
        end

        function output = get.PVY(obj)
            output = {obj.P, obj.V, obj.Y};
        end

        function output = get.PQ(obj)
            output = {obj.P, obj.Q};
        end

        function output = get.SH(obj)
            output = {obj.S, obj.H};
        end

        function output = get.SHX(obj)
            output = {obj.S, obj.H, obj.X};
        end

        function output = get.SHY(obj)
            output = {obj.S, obj.H, obj.Y};
        end

        function output = get.SP(obj)
            output = {obj.S, obj.P};
        end

        function output = get.SPX(obj)
            output = {obj.S, obj.P, obj.X};
        end

        function output = get.SPY(obj)
            output = {obj.S, obj.P, obj.Y};
        end

        function output = get.SPQ(obj)
            output = {obj.S, obj.P, obj.Q};
        end

        function output = get.ST(obj)
            output = {obj.S, obj.T};
        end

        function output = get.STX(obj)
            output = {obj.S, obj.T, obj.X};
        end

        function output = get.STY(obj)
            output = {obj.S, obj.T, obj.Y};
        end

        function output = get.SV(obj)
            output = {obj.S, obj.V};
        end

        function output = get.SVX(obj)
            output = {obj.S, obj.V, obj.X};
        end

        function output = get.SVY(obj)
            output = {obj.S, obj.V, obj.Y};
        end

        function output = get.SVQ(obj)
            output = {obj.S, obj.V, obj.Q};
        end

        function output = get.TD(obj)
            output = {obj.T, obj.D};
        end

        function output = get.TDX(obj)
            output = {obj.T, obj.D, obj.X};
        end

        function output = get.TDY(obj)
            output = {obj.T, obj.D, obj.Y};
        end

        function output = get.TDQ(obj)
            output = {obj.T, obj.D, obj.Q};
        end

        function output = get.TH(obj)
            output = {obj.T, obj.H};
        end

        function output = get.THX(obj)
            output = {obj.T, obj.H, obj.X};
        end

        function output = get.THY(obj)
            output = {obj.T, obj.H, obj.Y};
        end

        function output = get.TP(obj)
            output = {obj.T, obj.P};
        end

        function output = get.TPX(obj)
            output = {obj.T, obj.P, obj.X};
        end

        function output = get.TPY(obj)
            output = {obj.T, obj.P, obj.Y};
        end

        function output = get.TPQ(obj)
            output = {obj.T, obj.P, obj.Q};
        end

        function output = get.TQ(obj)
            output = {obj.T, obj.Q};
        end

        function output = get.TV(obj)
            output = {obj.T, obj.V};
        end

        function output = get.TVX(obj)
            output = {obj.T, obj.V, obj.X};
        end

        function output = get.TVY(obj)
            output = {obj.T, obj.V, obj.Y};
        end

        function output = get.UV(obj)
            output = {obj.U, obj.V};
        end

        function output = get.UVX(obj)
            output = {obj.U, obj.V, obj.X};
        end

        function output = get.UVY(obj)
            output = {obj.U, obj.V, obj.Y};
        end

        function output = get.UVQ(obj)
            output = {obj.U, obj.V, obj.Q};
        end

        function output = get.UP(obj)
            output = {obj.U, obj.P};
        end

        function output = get.UPX(obj)
            output = {obj.U, obj.P, obj.X};
        end

        function output = get.UPY(obj)
            output = {obj.U, obj.P, obj.Y};
        end

        function output = get.VH(obj)
            output = {obj.V, obj.H};
        end

        function output = get.VHX(obj)
            output = {obj.V, obj.H, obj.X};
        end

        function output = get.VHY(obj)
            output = {obj.V, obj.H, obj.Y};
        end

        %% Single-property setter methods

        function set.electricPotential(obj, phi)
            ctFunc('mThermo_setElectricPotential', obj.tpID, phi);
        end

        function set.basis(obj, b)

            if strcmp(b, 'mole') || strcmp(b, 'molar') ...
               || strcmp(b, 'Mole') || strcmp(b, 'Molar')
                obj.basis = 'molar';
            elseif strcmp(b, 'mass') || strcmp(b, 'Mass')
                obj.basis = 'mass';
            else
                error("Basis must be mass or molar.")
            end

        end

        function set.name(obj, str)
            ctFunc('mThermo_setName', obj.tpID, str);
        end

        function set.X(obj, xx)
            if isa(xx, 'double')
                ctFunc('mThermo_setMoleFractions', obj.tpID, xx);
            elseif isa(xx, 'char')
                ctFunc('mThermo_setMoleFractionsByName', obj.tpID, xx);
            else
                error('Invalid input.')
            end
        end

        function set.Y(obj, yy)
            if isa(yy, 'double')
                ctFunc('mThermo_setMassFractions', obj.tpID, yy);
            elseif isa(yy, 'char')
                ctFunc('mThermo_setMassFractionsByName', obj.tpID, yy);
            else
                error('Invalid input.')
            end

        end

        %% Multi-property setter methods

        function set.DP(obj, input)
            d = input{1};
            p = input{2};
            if strcmp(obj.basis, 'molar')
                d = d*obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_DP', obj.tpID, d, p);
        end

        function set.DPX(obj, input)
            obj.X = input{3};
            obj.DP = input(1:2);
        end

        function set.DPY(obj, input)
            obj.Y = input{3};
            obj.DP = input(1:2);
        end

        function set.HP(obj, input)
            h = input{1};
            p = input{2};
            if strcmp(obj.basis, 'molar')
                h = h/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_HP', obj.tpID, h, p);
        end

        function set.HPX(obj, input)
            obj.X = input{3};
            obj.HP = input(1:2);
        end

        function set.HPY(obj, input)
            obj.Y = input{3};
            obj.HP = input(1:2);
        end

        function set.PV(obj, input)
            p = input{1};
            v = input{2};
            if strcmp(obj.basis, 'molar')
                v = v/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_PV', obj.tpID, p, v);
        end

        function set.PVX(obj, input)
            obj.X = input{3};
            obj.PV = input(1:2);
        end

        function set.PVY(obj, input)
            obj.Y = input{3};
            obj.PV = input(1:2);
        end

        function set.PQ(obj, input)
            p = input{1};
            q = input{2};
            ctFunc('mThermo_setState_Psat', obj.tpID, p, q);
        end

        function set.SH(obj, input)
            s = input{1};
            h = input{2};
            if strcmp(obj.basis, 'molar')
                s = s/obj.meanMolecularWeight;
                h = h/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_SH', obj.tpID, s, h);
        end

        function set.SHX(obj, input)
            obj.X = input{3};
            obj.SH = input(1:2);
        end

        function set.SHY(obj, input)
            obj.Y = input{3};
            obj.SH = input(1:2);
        end

        function set.SP(obj, input)
            s = input{1};
            p = input{2};
            if strcmp(obj.basis, 'molar')
                s = s/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_SP', obj.tpID, s, p);
        end

        function set.SPX(obj, input)
            obj.X = input{3};
            obj.SP = input(1:2);
        end

        function set.SPY(obj, input)
            obj.Y = input{3};
            obj.SP = input(1:2);
        end

        function set.ST(obj, input)
            s = input{1};
            t = input{2};
            if strcmp(obj.basis, 'molar')
                s = s/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_ST', obj.tpID, s, t);
        end

        function set.STX(obj, input)
            obj.X = input{3};
            obj.ST = input(1:2);
        end

        function set.STY(obj, input)
            obj.Y = input{3};
            obj.ST = input(1:2);
        end

        function set.SV(obj, input)
            s = input{1};
            v = input{2};
            if strcmp(obj.basis, 'molar')
                s = s/obj.meanMolecularWeight;
                v = v/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_SV', obj.tpID, s, v);
        end

        function set.SVX(obj, input)
            obj.X = input{3};
            obj.SV = input(1:2);
        end

        function set.SVY(obj, input)
            obj.Y = input{3};
            obj.SV = input(1:2);
        end

        function set.TD(obj, input)
            t = input{1};
            d = input{2};
            if strcmp(obj.basis, 'molar')
                d = d*obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_TD', obj.tpID, t, d);
        end

        function set.TDX(obj, input)
            obj.X = input{3};
            obj.TD = input(1:2);
        end

        function set.TDY(obj, input)
            obj.Y = input{3};
            obj.TD = input(1:2);
        end

        function set.TH(obj, input)
            t = input{1};
            h = input{2};
            if strcmp(obj.basis, 'molar')
                h = h/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_TH', obj.tpID, t, h);
        end

        function set.THX(obj, input)
            obj.X = input{3};
            obj.TH = input(1:2);
        end

        function set.THY(obj, input)
            obj.Y = input{3};
            obj.TH = input(1:2);
        end

        function set.TP(obj, input)
            t = input{1};
            p = input{2};
            ctFunc('mThermo_setState_TP', obj.tpID, t, p);
        end

        function set.TPX(obj, input)
            t = input{1};
            p = input{2};
            x = input{3};
            if isa(x, 'char') || isa(x, 'string')
                ctFunc('mThermo_setState_TPX_byName', obj.tpID, t, p, x);
            elseif isa(x, 'double')
                ctFunc('mThermo_setState_TPX', obj.tpID, t, p, x);
            else
                error('Invalid input for mole fractions.')
            end
        end

        function set.TPY(obj, input)
            t = input{1};
            p = input{2};
            y = input{3};
            if isa(y, 'char') || isa(y, 'string')
                ctFunc('mThermo_setState_TPY_byName', obj.tpID, t, p, y);
            elseif isa(y, 'double')
                ctFunc('mThermo_setState_TPY', obj.tpID, t, p, y);
            else
                error('Invalid input for mass fractions.')
            end
        end

        function set.TQ(obj, input)
            t = input{1};
            q = input{2};
            ctFunc('mThermo_setState_Tsat', obj.tpID, t, q);
        end

        function set.TV(obj, input)
            t = input{1};
            v = input{2};
            if strcmp(obj.basis, 'molar')
                v = v/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_TV', obj.tpID, t, v);
        end

        function set.TVX(obj, input)
            obj.X = input{3};
            obj.TV = input(1:2);
        end

        function set.TVY(obj, input)
            obj.Y = input{3};
            obj.TV = input(1:2);
        end

        function set.UP(obj, input)
            u = input{1};
            p = input{2};
            if strcmp(obj.basis, 'molar')
                u = u/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_UP', obj.tpID, u, p);
        end

        function set.UPX(obj, input)
            obj.X = input{3};
            obj.UP = input(1:2);
        end

        function set.UPY(obj, input)
            obj.Y = input{3};
            obj.UP = input(1:2);
        end

        function set.UV(obj, input)
            u = input{1};
            v = input{2};
            if strcmp(obj.basis, 'molar')
                u = u/obj.meanMolecularWeight;
                v = v/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_UV', obj.tpID, u, v);
        end

        function set.UVX(obj, input)
            obj.X = input{3};
            obj.UV = input(1:2);
        end

        function set.UVY(obj, input)
            obj.Y = input{3};
            obj.UV = input(1:2);
        end

        function set.VH(obj, input)
            v = input{1};
            h = input{2};
            if strcmp(obj.basis, 'molar')
                v = v/obj.meanMolecularWeight;
                h = h/obj.meanMolecularWeight;
            end
            ctFunc('mThermo_setState_VH', obj.tpID, v, h);
        end

        function set.VHX(obj, input)
            obj.X = input{3};
            obj.VH = input(1:2);
        end

        function set.VHY(obj, input)
            obj.Y = input{3};
            obj.VH = input(1:2);
        end

        function setEquivalenceRatio(obj, phi, fuelComp, oxComp)
            % Set the mixture composition according to the equivalence ratio.
            ctFunc('mThermo_setEquivalenceRatio', obj.tpID, phi, fuelComp, oxComp);
        end
    end

end

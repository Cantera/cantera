classdef ThermoPhase < handle
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
    % :return:
    %     Instance of class :mat:class:`ThermoPhase`.

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
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase`.
        % :param b:
        %     String. Can be 'mole'/'molar'/'Molar'/'Mole' or 'mass'/'Mass'.
        basis

        name % Name of the phase.

        electricPotential % Electric potential. Units: V.

    end

    properties (SetAccess = protected)

        tpID % ID of the ThermoPhase object.

    end

    properties (SetAccess = immutable)

        atomicMasses % Atomic masses of the elements. Unit: kg/kmol.

        charges % Species charges. Unit: elem. charge.

        % Mean molecular weight ::
        %
        %     >> mmw = tp.meanMolecularWeight
        %
        % The mean molecular weight is the mole-fraction-weighted sum of the
        % molar masses of the individual species in the phase.
        %
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase` (or another
        %     object that derives from ThermoPhase).
        % :return:
        %     Scalar double mean molecular weight. Units: kg/kmol.
        meanMolecularWeight

        molarDensity % Molar basis density. Units: kmol/m^3.

        molecularWeights % Molecular weights of the species. Units: kg/kmol.

        nElements % Number of elements in the phase.

        nSpecies % Number of species in the phase.

        T % Temperature. Units: K.

        P % Pressure. Units: Pa.

        D % Density depending on the basis. Units: kmol/m^3 (molar) kg/m^3 (mass)

        H % Enthalpy depending on the basis. Units: J/kmol (molar) J/kg (mass).

        S % Entropy depending on the basis. Units: J/kmol-K (molar) J/kg-K (mass).

        U % Internal energy depending on the basis. Units: J/kmol (molar) J/kg (mass).

        G % Gibbs free energy depending on the basis. Units: J/kmol (molar) J/kg (mass).

        Q % Vapor fraction of the phase.

        % Basis-dependent specific heat at constant volume.
        % Units: J/kg-K (mass basis) or J/kmol-K (molar basis).
        cv

        % Basis-dependent specific heat at constant pressure.
        % Units: J/kg-K (mass basis) or J/kmol-K (molar basis).
        cp

        % Chemical potentials of the species. Units: J/kmol.
        chemicalPotentials

        % Electrochemical potentials of the species. Units: J/kmol.
        electrochemicalPotentials

        % Partial molar enthalpies for the species in the mixture. Units: J/kmol.
        partialMolarEnthalpies

        % Partial molar entropies for the species in the mixture. Units: J/kmol/K.
        partialMolarEntropies

        % Partial molar internal energies for the species in the mixture. Units: J/kmol.
        partialMolarIntEnergies

        % Partial molar heat capacities for the species in the mixture. Units: J/kmol/K.
        partialMolarCp

        % Partial molar volumes for the species in the mixture. Units: m^3/kmol.
        partialMolarVolumes

        critDensity % Critical density. Units: kg/m^3.

        critPressure % Critical pressure. Units: Pa.

        critTemperature % Critical temperature. Units: K.

        eosType % Type of equation of state.

        isIdealGas % A flag indicating whether the phase is an ideal gas.

        isothermalCompressibility % Isothermal compressibility. Units: 1/Pa.

        % Maximum temperature of the parameter fits ::
        %
        %     >> v = tp.maxTemp
        %
        % The parameterizations used to represent the temperature-dependent
        % species thermodynamic properties are generally only valid in some
        % finite temperature range, which may be different for each species
        % in the phase.
        %
        % See also: :mat:class:`minTemp`
        %
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase` (or another
        %     object that derives from ThermoPhase).
        % :return:
        %     Vector of maximum temperatures of all species.
        maxTemp

        % Minimum temperature of the parameter fits ::
        %
        %     >> v = tp.minTemp
        %
        % The parameterizations used to represent the temperature-dependent
        % species thermodynamic properties are generally only valid in some
        % finite temperature range, which may be different for each species
        % in the phase.
        %
        % See also: :mat:class:`maxTemp`
        %
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase` (or another
        %     object that derives from ThermoPhase).
        % :return:
        %     Vector of minimum temperatures of all species.
        minTemp

        refPressure % Reference pressure for standard-state. Units: Pa.

        % Generate a report describing the thermodynamic state of this phase.
        % To print the report to the terminal, simply call the phase object.
        % The following two statements are equivalent ::
        %
        %     >> phase
        %     >> disp(phase.report)
        report

        % Saturation pressure at current temperature ::
        %
        %     >> p = tp.satPressure.
        %
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase` (or another
        %     object that derives from ThermoPhase).
        % :return:
        %     Saturation pressure for temperature T. Units: Pa.
        satPressure

        % Saturation temperature at current pressure ::
        %
        %     >> t = tp.satTemperature.
        %
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase` (or another
        %     object that derives from ThermoPhase).
        % :return:
        %     Saturation temperature for pressure p. Units: K.
        satTemperature

        % Speed of sound ::
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
        %
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase` (or another
        %     class derived from ThermoPhase).
        % :return:
        %     The speed of sound. Units: m/s.
        soundSpeed

        % All species names ::
        %
        %     >> n = tp.speciesNames
        %
        % :param tp:
        %     Instance of class :mat:class:`ThermoPhase` (or another
        %     class derived from ThermoPhase).
        % :return:
        %     Cell array of strings of all of the species names.
        speciesNames

        thermalExpansionCoeff % Thermal expansion coefficient. Units: 1/K.

    end

    properties (Dependent = true)

        V % Basis-dependent specific volume. Units: m^3/kmol (molar) m^3/kg (mass).

        % Get/Set density [kg/m^3 or kmol/m^3] and pressure [Pa].
        DP

        % Get/Set density [kg/m^3 or kmol/m^3], pressure [Pa], and mole fractions.
        DPX

        % Get/Set density [kg/m^3 or kmol/m^3], pressure [Pa], and mass fractions.
        DPY

        % Get density [kg/m^3 or kmol/m^3], pressure [Pa], and vapor fraction.
        DPQ

        % Get/Set enthalpy [J/kg or J/kmol] and pressure [Pa].
        HP

        % Get/Set enthalpy [J/kg or J/kmol], pressure [Pa], and mole fractions.
        HPX

        % Get/Set enthalpy [J/kg or J/kmol], pressure [Pa], and mass fractions.
        HPY

        % Get enthalpy [J/kg or J/kmol], pressure [Pa], and vapor fraction.
        HPQ

        % Get/Set pressure [Pa] and specific volume [m^3/kg or m^3/kmol].
        PV

        % Get/Set pressure [Pa], specific volume [m^3/kg or m^3/kmol],
        % and mole fractions.
        PVX

        % Get/Set pressure [Pa], specific volume [m^3/kg or m^3/kmol],
        % and mass fractions.
        PVY

        % Get/Set pressure [Pa] and vapor fraction of a two-phase state.
        PQ

        % Get/Set entropy [J/kg-K or J/kmol-K] and enthalpy [J/kg or J/kmol].
        SH

        % Get/Set entropy [J/kg-K or J/kmol-K], enthalpy [J/kg or J/kmol],
        % and mole fractions.
        SHX

        % Get/Set entropy [J/kg-K or J/kmol-K], enthalpy [J/kg or J/kmol],
        % and mass fractions.
        SHY

        % Get/Set entropy [J/kg-K or J/kmol-K] and pressure [Pa].
        SP

        % Get/Set entropy [J/kg-K or J/kmol-K], pressure [Pa], and mole fractions.
        SPX

        % Get/Set entropy [J/kg-K or J/kmol-K], pressure [Pa], and mass fractions.
        SPY

        % Get entropy [J/kg-K or J/kmol-K], pressure [Pa], and vapor fraction.
        SPQ

        % Get/Set entropy [J/kg-K or J/kmol-K] and temperature [K].
        ST

        % Get/Set entropy [J/kg-K or J/kmol-K], temperature [K], and mole fractions.
        STX

        % Get/Set entropy [J/kg-K or J/kmol-K], temperature [K], and mass fractions.
        STY

        % Get/Set entropy [J/kg-K or J/kmol-K] and specific volume [m^3/kg or m^3/kmol].
        SV

        % Get/Set entropy [J/kg-K or J/kmol-K], specific volume [m^3/kg or m^3/kmol],
        % and mole fractions.
        SVX

        % Get/Set entropy [J/kg-K or J/kmol-K], specific volume [m^3/kg or m^3/kmol],
        % and mass fractions.
        SVY

        % Get/Set entropy [J/kg-K or J/kmol-K], specific volume [m^3/kg or m^3/kmol],
        % and vapor fraction.
        SVQ

        % Get/Set temperature [K] and density [kg/m^3 or kmol/m^3].
        TD

        % Get/Set temperature [K], density [kg/m^3 or kmol/m^3], and mole fractions.
        TDX

        % Get/Set temperature [K], density [kg/m^3 or kmol/m^3], and mass fractions.
        TDY

        % Get temperature [K], density [kg/m^3 or kmol/m^3], and vapor fraction.
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

        % Get/Set temperature [K] and specific volume [m^3/kg or m^3/kmol].
        TV

        % Get/Set temperature [K], specific volume [m^3/kg or m^3/kmol],
        % and mole fractions.
        TVX

        % Get/Set temperature [K], specific volume [m^3/kg or m^3/kmol],
        % and mass fractions.
        TVY

        % Get/Set internal energy [J/kg or J/kmol] and specific volume
        % [m^3/kg or m^3/kmol].
        UV

        % Get/Set internal energy [J/kg or J/kmol], specific volume
        % [m^3/kg or m^3/kmol], and mole fractions.
        UVX

        % Get/Set internal energy [J/kg or J/kmol], specific volume
        % [m^3/kg or m^3/kmol], and mass fractions.
        UVY

        % Get internal energy [J/kg or J/kmol], specific volume
        % [m^3/kg or m^3/kmol], and vapor fraction.
        UVQ

        % Get/Set internal energy [J/kg or J/kmol] and pressure [Pa].
        UP

        % Get/Set internal energy [J/kg or J/kmol], pressure [Pa],
        % and mole fractions.
        UPX

        % Get/Set internal energy [J/kg or J/kmol], pressure [Pa],
        % and mass fractions.
        UPY

        % Get/Set volume [m^3/kg or m^3/kmol] and enthalpy [J/kg or J/kmol].
        VH

        % Get/Set specific volume [m^3/kg or m^3/kmol], enthalpy [J/kg or J/kmol],
        % and mole fractions.
        VHX

        % Get/Set specific volume [m^3/kg or m^3/kmol], enthalpy [J/kg or J/kmol],
        % and mass fractions.
        VHY

    end

    methods
        %% ThermoPhase Class Constructor

        function tp = ThermoPhase(id)
            % Create a :mat:class:`ThermoPhase` object.
            ctIsLoaded;

            if ~isnumeric(id)
                error('Invalid argument: constructor requires integer solution ID.')
            end

            tp.tpID = ctFunc('soln_thermo', id);
            tp.basis = 'molar';
        end

        %% ThermoPhase Utility Methods

        function display(tp)
            disp(tp.report);
        end

        function tp = equilibrate(tp, xy, solver, rtol, maxsteps, maxiter, loglevel)
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
            %     Specifies the equilibrium solver to use. If solver = 0, a fast
            %     solver using the element potential method will be used. If
            %     solver = 1, a slower but more robust Gibbs minimization solver
            %     will be used. If solver >= 2, a version of the VCS algorithm will
            %     be used. If solver < 0 or is unspecified, the fast solver
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

            ctFunc('thermo_equilibrate', tp.tpID, xy, solver, rtol, ...
                    maxsteps, maxiter, loglevel);
        end

        %% ThermoPhase inquiry methods

        function k = elementIndex(tp, name)
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
                        k(i, j) = ctFunc('thermo_elementIndex', ...
                                        tp.tpID, name{i, j}) + 1;

                        if k(i, j) > 1e3
                            warning(['Element ', name{i, j}, ...
                                    ' does not exist in the phase']);
                            k(i, j) = -1;
                        end

                    end

                end

            elseif ischar(name)
                k = ctFunc('thermo_elementIndex', tp.tpID, name) + 1;

                if k > 1e3
                    warning(['Element ', name, ' does not exist in the phase']);
                    k = -1;
                end

            else
                error('name must be either a string or cell array of strings')
            end

        end

        function elMassFrac = elementalMassFraction(tp, element)
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
            % :param tp:
            %     Object representing the gas, instance of class :mat:class:`Solution`,
            %     and an ideal gas. The state of this object should be set to an
            %     estimate of the gas state before calling `elementalMassFraction`.
            % :param element:
            %     String representing the element name.
            % :return:
            %     Elemental mass fraction within a gas object.

            if nargin ~= 2
                error('elementalMassFraction expects two input arguments.');
            end

            if ~tp.isIdealGas
                error('Gas object must represent an ideal gas mixture.');
            end

            if ~ischar(element)
                error('Wrong type for element name: must be character.');
            end

            n = tp.nSpecies;
            spec = tp.speciesNames;
            eli = tp.elementIndex(element);
            M = tp.atomicMasses;
            Mel = M(eli);
            MW = tp.molecularWeights;
            % Initialize the element mass fraction as zero.
            elMassFrac = 0.0;
            % Perform summation of elemental mass fraction over all species.
            for i = 1:n
                natoms(i) = tp.nAtoms(spec{i}, element);
                yy(i) = tp.massFraction(spec{i});
                elMassFrac = elMassFrac + (natoms(i) * Mel * yy(i)) / MW(i);
            end

        end

        function n = nAtoms(tp, species, element)
            % Number of atoms of an element in a species ::
            %
            %   >> n = tp.nAtoms(k,m)
            %
            % :param tp:
            %     Instance of class :mat:class:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param k:
            %     String species name or integer species number
            % :param m:
            %     String element name or integer element number
            % :return:
            %     Number of atoms of element ``m`` in species ``k``.

            if nargin ~= 3
                error('Two input arguments required.')
            end

            if ischar(species)
                k = tp.speciesIndex(species);
            else
                k = species;
            end

            if ischar(element)
                m = tp.elementIndex(element);
            else
                m = element;
            end

            n = ctFunc('thermo_nAtoms', tp.tpID, k - 1, m - 1);

        end

        function k = speciesIndex(tp, name)
            % Index of a species given the name ::
            %
            %   >> k = tp.speciesIndex(name)
            %
            % The index is an integer assigned to each species in sequence as it
            % is read in from the input file.
            %
            % NOTE: In keeping with the conventions used by Matlab, this method
            % returns 1 for the first species, 2 for the second, etc. In
            % contrast, the corresponding method speciesIndex in the Cantera C++
            % and Python interfaces returns 0 for the first species, 1 for the
            % second one, etc. ::
            %
            %     >> ich4 = gas.speciesIndex('CH4');
            %     >> iho2 = gas.speciesIndex('HO2');
            %
            % :param tp:
            %     Instance of class :mat:class:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
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
                        k(i, j) = ctFunc('thermo_speciesIndex', ...
                                        tp.tpID, name{i, j}) + 1;

                        if k(i, j) > 1e6
                            warning(['Species ', name{i, j}, ...
                                    ' does not exist in the phase']);
                            k(i, j) = -1;
                        end

                    end

                end

            elseif ischar(name)
                k = ctFunc('thermo_speciesIndex', tp.tpID, name) + 1;

                if k > 1e6
                    warning(['Species ', name, ' does not exist in the phase.']);
                    k = -1;
                end

            else
                error('Name must be either a string or cell array of strings.')
            end

        end

        function nm = speciesName(tp, k)
            % Name of one or multiple species given the index ::
            %
            %     >> k = tp.speciesName(k)
            %
            % :param tp:
            %     Instance of class :mat:class:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
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
                    output = ctString('thermo_getSpeciesName', tp.tpID, ksp);
                    nm{i, j} = output;
                end

            end

        end

        function x = moleFraction(tp, species)
            % Get the mole fraction of one or a list of species ::
            %
            %     >> x = tp.moleFraction(species)
            %
            % :param tp:
            %     Instance of class :mat:class:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param species:
            %     String or cell array of strings of species whose mole
            %     fraction is desired
            % :return:
            %     Scalar or vector double mole fractions

            xarray = tp.X;

            if isa(species, 'char')
                k = tp.speciesIndex(species);

                if k > 0
                    x = xarray(k);
                else
                    error("species not found.");
                end

            elseif isa(species, 'cell')
                n = length(species);
                x = zeros(1, n);

                for j = 1:n
                    k = tp.speciesIndex(species{j});

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

        function y = massFraction(tp, species)
            % Get the mass fraction of one or a list of species ::
            %
            %     >> y = tp.massFraction(species)
            %
            % :param tp:
            %     Instance of class :mat:class:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param species:
            %     String or cell array of strings of species whose mass
            %     fraction is desired
            % :return:
            %     Scalar or vector double mass fractions

            yy = tp.Y;

            if isa(species, 'char')
                k = tp.speciesIndex(species);

                if k > 0
                    y = yy(k);
                else
                    error("Error: species not found.")
                end

            elseif isa(species, 'cell')
                n = length(species);
                y = zeros(1, n);

                for j = 1:n
                    k = tp.speciesIndex(species{j});

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

        %% Single-property getter methods

        function amu = get.atomicMasses(tp)
            nel = tp.nElements;
            aa = zeros(1, nel);
            pt = libpointer('doublePtr', aa);
            ctFunc('thermo_getAtomicWeights', tp.tpID, nel, pt);
            amu = pt.Value;
        end

        function e = get.charges(tp)
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            ctFunc('thermo_getCharges', tp.tpID, nsp, pt);
            e = pt.Value;
        end

        function c = get.cv(tp)
            if strcmp(tp.basis, 'molar')
                c = ctFunc('thermo_cv_mole', tp.tpID);
            else
                c = ctFunc('thermo_cv_mass', tp.tpID);
            end
        end

        function c = get.cp(tp)
            if strcmp(tp.basis, 'molar')
                c = ctFunc('thermo_cp_mole', tp.tpID);
            else
                c = ctFunc('thermo_cp_mass', tp.tpID);
            end
        end

        function d = get.critDensity(tp)
            d = ctFunc('thermo_critDensity', tp.tpID);
        end

        function p = get.critPressure(tp)
            p = ctFunc('thermo_critPressure', tp.tpID);
        end

        function t = get.critTemperature(tp)
            t = ctFunc('thermo_critTemperature', tp.tpID);
        end

        function v = get.electricPotential(tp)
            v = ctFunc('thermo_electricPotential', tp.tpID);
        end

        function e = get.eosType(tp)
            e = ctString('thermo_getEosType', tp.tpID);
        end

        function v = get.isIdealGas(tp)
            v = strcmp(tp.eosType, 'ideal-gas');
        end

        function b = get.isothermalCompressibility(tp)
            b = ctFunc('thermo_isothermalCompressibility', tp.tpID);
        end

        function t = get.maxTemp(tp)
            t = ctFunc('thermo_maxTemp', tp.tpID, -1);
        end

        function t = get.minTemp(tp)
            t = ctFunc('thermo_minTemp', tp.tpID, -1);
        end

        function mmw = get.meanMolecularWeight(tp)
            mmw = ctFunc('thermo_meanMolecularWeight', tp.tpID);
        end

        function density = get.molarDensity(tp)
            density = ctFunc('thermo_molarDensity', tp.tpID);
        end

        function mw = get.molecularWeights(tp)
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            ctFunc('thermo_getMolecularWeights', tp.tpID, nsp, pt);
            mw = pt.Value;
        end

        function nel = get.nElements(tp)
            nel = ctFunc('thermo_nElements', tp.tpID);
        end

        function nsp = get.nSpecies(tp)
            nsp = ctFunc('thermo_nSpecies', tp.tpID);
        end

        function p = get.refPressure(tp)
            p = ctFunc('thermo_refPressure', tp.tpID);
        end

        function p = get.satPressure(tp)
            p = ctFunc('thermo_satPressure', tp.tpID, tp.T);
        end

        function t = get.satTemperature(tp)
            t = ctFunc('thermo_satTemperature', tp.tpID, tp.P);
        end

        function c = get.soundSpeed(tp)

            if tp.isIdealGas
                tp.basis = 'mass';
                gamma = tp.cp / tp.cv;
                wtm = tp.meanMolecularWeight;
                r = 8314.4621 / wtm;
                c = sqrt(gamma * r * tp.T);
            else
                rho0 = tp.D;
                p0 = tp.P;
                s0 = tp.S;
                rho1 = 1.001 * rho0;
                tp.SV = {s0, 1/rho0};
                p1 = tp.P;
                dpdrho_s = (p1 - p0) / (rho1 - rho0);
                c = sqrt(dpdrho_s);
            end

        end

        function s = get.name(tp)
            s = ctString('thermo_getName', tp.tpID);
        end

        function str = get.report(tp)
            buflen = 0 - calllib(ctLib, 'thermo_report', tp.tpID, 0, '', 1);
            aa = char(ones(1, buflen));
            ptr = libpointer('cstring', aa);
            [iok, bb] = calllib(ctLib, 'thermo_report', tp.tpID, buflen, ptr, 1);

            if iok < 0
                error(ctGetErr);
            end

            str = bb;
        end

        function n = get.speciesNames(tp)
            n = tp.speciesName(1:tp.nSpecies);
        end

        function a = get.thermalExpansionCoeff(tp)
            a = ctFunc('thermo_thermalExpansionCoeff', tp.tpID);
        end

        function temperature = get.T(tp)
            temperature = ctFunc('thermo_temperature', tp.tpID);
        end

        function pressure = get.P(tp)
            pressure = ctFunc('thermo_pressure', tp.tpID);
        end

        function v = get.Q(tp)
            v = ctFunc('thermo_vaporFraction', tp.tpID);
        end

        function density = get.D(tp)
            if strcmp(tp.basis, 'mass')
                density = ctFunc('thermo_density', tp.tpID);
            else
                density = tp.molarDensity;
            end
        end

        function volume = get.V(tp)
            volume = 1 / tp.D;
        end

        function moleFractions = get.X(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_getMoleFractions', tp.tpID, nsp, pt);
            moleFractions = pt.Value;
        end

        function massFractions = get.Y(tp)
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            ctFunc('thermo_getMassFractions', tp.tpID, nsp, pt);
            massFractions = pt.Value;
        end

        function enthalpy = get.H(tp)
            if strcmp(tp.basis, 'molar')
                enthalpy = ctFunc('thermo_enthalpy_mole', tp.tpID);
            else
                enthalpy = ctFunc('thermo_enthalpy_mass', tp.tpID);
            end
        end

        function mu = get.chemicalPotentials(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_chemPotentials', tp.tpID, nsp, pt);
            mu = pt.Value;
        end

        function emu = get.electrochemicalPotentials(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_electrochemPotentials', tp.tpID, nsp, pt);
            emu = pt.Value;
        end

        function enthalpies = get.partialMolarEnthalpies(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_getPartialMolarEnthalpies', tp.tpID, nsp, pt);
            enthalpies = pt.Value;
        end

        function entropies = get.partialMolarEntropies(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_getPartialMolarEntropies', tp.tpID, nsp, pt);
            entropies = pt.Value;
        end

        function intEnergies = get.partialMolarIntEnergies(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_getPartialMolarIntEnergies', tp.tpID, nsp, pt);
            intEnergies = pt.Value;
        end

        function cps = get.partialMolarCp(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_getPartialMolarCp', tp.tpID, nsp, pt);
            cps = pt.Value;
        end

        function volumes = get.partialMolarVolumes(tp)
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('thermo_getPartialMolarVolumes', tp.tpID, nsp, pt);
            volumes = pt.Value;
        end

        function entropy = get.S(tp)
            if strcmp(tp.basis, 'molar')
                entropy = ctFunc('thermo_entropy_mole', tp.tpID);
            else
                entropy = ctFunc('thermo_entropy_mass', tp.tpID);
            end
        end

        function intEnergy = get.U(tp)
            if strcmp(tp.basis, 'molar')
                intEnergy = ctFunc('thermo_intEnergy_mole', tp.tpID);
            else
                intEnergy = ctFunc('thermo_intEnergy_mass', tp.tpID);
            end
        end

        function gibbs = get.G(tp)
            if strcmp(tp.basis, 'molar')
                gibbs = ctFunc('thermo_gibbs_mole', tp.tpID);
            else
                gibbs = ctFunc('thermo_gibbs_mass', tp.tpID);
            end
        end

        %% Multi-property getter methods

        function output = get.DP(tp)
            output = {tp.D, tp.P};
        end

        function output = get.DPX(tp)
            output = {tp.D, tp.P, tp.X};
        end

        function output = get.DPY(tp)
            output = {tp.D, tp.P, tp.Y};
        end

        function output = get.DPQ(tp)
            output = {tp.D, tp.P, tp.Q};
        end

        function output = get.HP(tp)
            output = {tp.H, tp.P};
        end

        function output = get.HPX(tp)
            output = {tp.H, tp.P, tp.X};
        end

        function output = get.HPY(tp)
            output = {tp.H, tp.P, tp.Y};
        end

        function output = get.HPQ(tp)
            output = {tp.H, tp.P, tp.Q};
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

        function output = get.PQ(tp)
            output = {tp.P, tp.Q};
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

        function output = get.SPQ(tp)
            output = {tp.S, tp.P, tp.Q};
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

        function output = get.SVQ(tp)
            output = {tp.S, tp.V, tp.Q};
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

        function output = get.TDQ(tp)
            output = {tp.T, tp.D, tp.Q};
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
            output = {tp.T, tp.P, tp.X};
        end

        function output = get.TPY(tp)
            output = {tp.T, tp.P, tp.Y};
        end

        function output = get.TPQ(tp)
            output = {tp.T, tp.P, tp.Q};
        end

        function output = get.TQ(tp)
            output = {tp.T, tp.Q};
        end

        function output = get.TV(tp)
            output = {tp.T, tp.V};
        end

        function output = get.TVX(tp)
            output = {tp.T, tp.V, tp.X};
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

        function output = get.UVQ(tp)
            output = {tp.U, tp.V, tp.Q};
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

        %% Single-property setter methods

        function tp = set.electricPotential(tp, phi)
            ctFunc('thermo_setElectricPotential', tp.tpID, phi);
        end

        function set.basis(tp, b)

            if strcmp(b, 'mole') || strcmp(b, 'molar') ...
                || strcmp(b, 'Mole') || strcmp(b, 'Molar')
                tp.basis = 'molar';
            elseif strcmp(b, 'mass') || strcmp(b, 'Mass')
                tp.basis = 'mass';
            else
                error("Basis must be mass or molar.")
            end

        end

        function set.name(tp, str)
            ctFunc('thermo_setName', tp.tpID, str);
        end

        function set.X(tp, xx)
            tol = 1e-9;

            if isempty(xx)
                error('Array cannot be empty');
            end

            if isa(xx, 'double')
                nsp = tp.nSpecies;
                if length(xx) ~= nsp
                    error('Length of array must be equal to number of species.')
                end

                if length(xx) ~= nsp
                    error('Length of array must be equal to number of species.')
                end

                if abs(sum(xx) - 1) <= tol
                    norm = 0;
                else
                    norm = 1;
                end

                ctFunc('thermo_setMoleFractions', tp.tpID, nsp, xx, norm);
            elseif isa(xx, 'char')
                ctFunc('thermo_setMoleFractionsByName', tp.tpID, xx);
            else
                error('Invalid input.')
            end

        end

        function set.Y(tp, yy)
            tol = 1e-9;

            if isempty(yy)
                error('Array cannot be empty');
            end

            if isa(yy, 'double')
                nsp = tp.nSpecies;

                if length(yy) ~= nsp
                    error('Length of array must be equal to number of species.')
                end

                if abs(sum(yy) -1) <= tol
                    norm = 0;
                else
                    norm = 1;
                end

                ctFunc('thermo_setMassFractions', tp.tpID, nsp, yy, norm);
            elseif isa(yy, 'char')
                ctFunc('thermo_setMassFractionsByName', tp.tpID, yy);
            else
                error('Invalid input.')
            end

        end

        %% Multi-property setter methods

        function set.DP(tp, input)
            d = input{1};
            p = input{2};
            if strcmp(tp.basis, 'molar')
                d = d*tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_DP', tp.tpID, [d, p]);
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
            h = input{1};
            p = input{2};
            if strcmp(tp.basis, 'molar')
                h = h/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_HP', tp.tpID, [h, p]);
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
            p = input{1};
            v = input{2};
            if strcmp(tp.basis, 'molar')
                v = v/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_PV', tp.tpID, [p, v]);
        end

        function set.PVX(tp, input)
            tp.X = input{3};
            tp.PV = input(1:2);
        end

        function set.PVY(tp, input)
            tp.Y = input{3};
            tp.PV = input(1:2);
        end

        function set.PQ(tp, input)
            p = input{1};
            q = input{2};
            ctFunc('thermo_setState_Psat', tp.tpID, p, q);
        end

        function set.SH(tp, input)
            s = input{1};
            h = input{2};
            if strcmp(tp.basis, 'molar')
                s = s/tp.meanMolecularWeight;
                h = h/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_SH', tp.tpID, [s, h]);
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
            s = input{1};
            p = input{2};
            if strcmp(tp.basis, 'molar')
                s = s/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_SP', tp.tpID, [s, p]);
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
            s = input{1};
            t = input{2};
            if strcmp(tp.basis, 'molar')
                s = s/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_ST', tp.tpID, [s, t]);
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
            s = input{1};
            v = input{2};
            if strcmp(tp.basis, 'molar')
                s = s/tp.meanMolecularWeight;
                v = v/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_SV', tp.tpID, [s, v]);
        end

        function set.SVX(tp, input)
            tp.X = input{3};
            tp.SV = input(1:2);
        end

        function set.SVY(tp, input)
            tp.Y = input{3};
            tp.SV = input(1:2);
        end

        function set.TD(tp, input)
            t = input{1};
            d = input{2};
            if strcmp(tp.basis, 'molar')
                d = d*tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_TD', tp.tpID, [t, d]);
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
            t = input{1};
            h = input{2};
            if strcmp(tp.basis, 'molar')
                h = h/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_TH', tp.tpID, [t, h]);
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
            ctFunc('thermo_set_TP', tp.tpID, [t, p]);
        end

        function set.TPX(tp, input)
            tp.X = input{3};
            tp.TP = input(1:2);
        end

        function set.TPY(tp, input)
            tp.Y = input{3};
            tp.TP = input(1:2);
        end

        function tp = set.TQ(tp, input)
            t = input{1};
            q = input{2};
            ctFunc('thermo_setState_Tsat', tp.tpID, t, q);
        end

        function set.TV(tp, input)
            t = input{1};
            v = input{2};
            if strcmp(tp.basis, 'molar')
                v = v/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_TV', tp.tpID, [t, v]);
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
            u = input{1};
            p = input{2};
            if strcmp(tp.basis, 'molar')
                u = u/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_UP', tp.tpID, [u, p]);
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
            u = input{1};
            v = input{2};
            if strcmp(tp.basis, 'molar')
                u = u/tp.meanMolecularWeight;
                v = v/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_UV', tp.tpID, [u, v]);
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
            v = input{1};
            h = input{2};
            if strcmp(tp.basis, 'molar')
                v = v/tp.meanMolecularWeight;
                h = h/tp.meanMolecularWeight;
            end
            ctFunc('thermo_set_VH', tp.tpID, [v, h]);
        end

        function set.VHX(tp, input)
            tp.X = input{3};
            tp.VH = input(1:2);
        end

        function set.VHY(tp, input)
            tp.Y = input{3};
            tp.VH = input(1:2);
        end

    end

end

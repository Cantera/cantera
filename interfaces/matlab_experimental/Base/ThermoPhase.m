classdef ThermoPhase < handle

    properties
        tpID
        T % temperature
        P % pressure
        D % density
        X % mole fractions
        Y % mass fractions
        H % enthalpy
        S % entropy
        U % internal energy
        G % Gibbs free energy
        basis
    end

    properties (Dependent)
        V % specific volume
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
            % THERMOPHASE  ThermoPhase class constructor.
            % t = ThermoPhase(src, id)
            % :param src:
            %     Input string of YAML, CTI, or XML file name.
            % :param id:
            %     ID of the phase to import as specified in the input file. (optional)
            % :return:
            %     Instance of class :mat:func:`ThermoPhase`
            %
            checklib;
            if nargin > 2
                error('ThermoPhase expects 1 or 2 input arguments.');
            end
            if nargin == 1
                id = '-';
            end
            tp.tpID = calllib(ct, 'thermo_newFromFile', src, id);
            tp.basis = 'molar';
        end

        %% Utility methods

        function display(tp, threshold)
            % Display thermo properties

            if nargin < 2 || ~isnumeric(threshold)
                threshold = 1e-14;
            end
            calllib(ct, 'thermo_print', tp.tpID, 1, threshold);
        end

        function tpClear(tp)
            % CLEAR  Delete the kernel object.
            % tp.tpClear
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            %
            calllib(ct, 'thermo_del', tp.tpID);
        end

        function tp = set.basis(tp, b)
            % BASIS    Determines whether intensive thermodynamic properties are
            % treated on a mass (per kg) or molar (per kmol) basis. This
            % affects the values returned by the properties H, U, S, G, V,
            % Density, Cv, and Cp, as well as the values used with the
            % state-setting properties such as HPX and UV.
            % tp.Basis
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase`.
            % :param b:
            %     String. Can be 'mole'/'molar'/'Molar'/Mole' or 'mass'/'Mass'.

            if strcmp(b, 'mole') || strcmp(b, 'molar') ...
               || strcmp(b, 'Mole') || strcmp(b, 'Molar')
                tp.basis = 'molar';
            elseif strcmp(b, 'mass') || strcmp(b, 'Mass')
                tp.basis = 'mass';
            else error("Basis must be mass or molar")
            end
        end

        function tp = equilibrate(tp, xy, solver, rtol, maxsteps, ...
                                  maxiter, loglevel)
            % EQUILIBRATE  Set the phase to a state of chemical equilibrium.
            % tp.equilibrate(xy, solver, rtol, maxsteps, maxiter, loglevel)
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
            %
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
            calllib(ct, 'thermo_equilibrate', tp.tpID, xy, solver, ...
                    rtol, maxsteps, maxiter, loglevel);
        end

        %% PhaseGet single methods

        function amu = atomicMasses(tp)
            % ATOMICMASSES  Get the atomic masses of the elements.
            % x = tp.atomicMasses
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase).
            % :return:
            %     Vector of element atomic masses. Units: kg/kmol
            %
            nel = tp.nElements;
            aa = zeros(1, nel);
            pt = libpointer('doublePtr', aa);
            calllib(ct, 'thermo_getAtomicWeights', ...
                    tp.tpID, nel, pt);
            amu = pt.Value;
        end

        function e = charges(tp)
            % CHARGES  Get the array of species charges
            % x = tp.charges
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Vector of species charges. Units: elem. charge
            %
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(ct, 'thermo_getCharges', ...
                    tp.tpID, nsp, pt);
            e = pt.Value;
        end

        function k = elementIndex(tp, name)
            % ELEMENTINDEX  Get the index of an element given its name.
            % k = tp.elementIndex(name)
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
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param name:
            %     String or cell array of strings of elements to look up
            % :return:
            %     Integer or vector of integers of element indices
            %
            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);
                for i = 1:m
                    for j = 1:n
                        k(i, j) = calllib(ct, 'thermo_elementIndex', ...
                                          tp.tpID, name{i, j}) + 1;
                        if k(i, j) > 1e3
                            warning(['Element ', name{i, j}, ...
                                   ' does not exist in the phase']);
                            k(i, j) = -1;
                        end
                    end
                end
            elseif ischar(name)
                k = calllib(ct, 'thermo_elementIndex', ...
                            tp.tpID, name) + 1;
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
            % ELEMENTALMASSFRACTION  Determine the elemental mass fraction in gas object.
            % elMassFrac = tp.elementalMassFraction(element)
            % :param tp:
            %     Object representing the gas, instance of class :mat:func:`Solution`,
            %     and an ideal gas. The state of this object should be set to an
            %     estimate of the gas state before calling elementalMassFraction.
            % :param element:
            %     String representing the element name.
            % :return:
            %     Elemental mass fraction within a gas object.
            %
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
            % MEANMOLECULARWEIGHT  Get the mean molecular weight.
            % mmw = tp.meanMolecularWeight
            % The mean molecular weight is the mole-fraction-weighted sum of the
            % molar masses of the individual species in the phase.
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Scalar double mean molecular weight. Units: kg/kmol
            %
            mmw = calllib(ct, 'thermo_meanMolecularWeight', tp.tpID);
        end

        function density = molarDensity(tp)
            % Get the molar basis density in kmol/m^3.

            density = calllib(ct, 'thermo_molarDensity', tp.tpID);
        end

        function mw = MolecularWeights(tp)
            % MOLECULARWEIGHTS  Get the molecular weights of the species.
            % mw = tp.molecularWeights
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Vector of species molecular weights. Units: kg/kmol
            %
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(ct, 'thermo_getMolecularWeights', ...
                    tp.tpID, nsp, pt);
            mw = pt.Value;
        end

        function n = nAtoms(tp, species, element)
            % NATOMS  Get the number of atoms of an element in a species.
            % n = tp.nAtoms(k,m)
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param k:
            %     String species name or integer species number
            % :param m:
            %     String element name or integer element number
            % :return:
            %     Number of atoms of element ``m`` in species ``k``.
            %
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
                n = calllib(ct, 'thermo_nAtoms', tp.tpID, k-1, m-1);
            else
                error('Two input arguments required.')
            end
        end

        function nel = nElements(tp)
            % NELEMENTS  Get the number of elements.
            % n = tp.nElements
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Number of elements in the phase.
            %
            nel = calllib(ct, 'thermo_nElements', tp.tpID);
        end

        function nsp = nSpecies(tp)
            % NSPECIES  Get the number of species.
            % n = tp.nSpecies
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Number of species in the phase.
            %
            nsp = calllib(ct, 'thermo_nSpecies', tp.tpID);
        end

        function k = speciesIndex(tp, name)
            % SPECIESINDEX  Get the index of a species given the name.
            % k = tp.speciesIndex(name)
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
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param name:
            %     If name is a single string, the return value will be a integer
            %     containing the corresponding index. If it is an cell array of
            %     strings, the output will be an array of the same shape
            %     containing the indices.
            % :return:
            %     Scalar or array of integers
            %
            if iscell(name)
                [m, n] = size(name);
                k = zeros(m, n);
                for i = 1:m
                    for j = 1:n
                        k(i, j) = calllib(ct, 'thermo_speciesIndex', ...
                                          tp.tpID, name{i, j}) + 1;
                        if k(i, j) > 1e3
                            warning(['Species ', name{i, j}, ...
                                   ' does not exist in the phase']);
                            k(i, j) = -1;
                        end
                    end
                end
            elseif ischar(name)
                k = calllib(ct, 'thermo_speciesIndex', ...
                            tp.tpID, name) + 1;
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
            % SPECIESNAME  Get the name of a species given the index.
            % nm = tp.speciesName(k)
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param k:
            %     Scalar or array of integer species numbers
            % :return:
            %     Cell array of strings
            %
            [m, n] = size(k);
            nm = cell(m, n);
            for i = 1:m
                for j = 1:n
                    ksp = k(i, j) - 1;
                    buflen = calllib(ct, 'thermo_getSpeciesName', ...
                                     tp.tpID, ksp, 0, '');
                    if buflen > 0
                        aa = char(zeros(1, buflen));
                        [~, aa] = calllib(ct, 'thermo_getSpeciesName', ...
                                          tp.tpID, ksp, buflen, aa);
                        nm{i, j} = aa;
                    end
                end
            end
        end

        function n = speciesNames(tp)
            % SPECIESNAMES  Get the species names.
            % n = tp.speciesNames
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :return:
            %     Cell array of strings of all of the species names
            %
            n = tp.speciesName(1:tp.nSpecies);
        end

        function temperature = get.T(tp)
            % GET.T  Get the temperature.
            % temperature = tp.T
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :return:
            %     Temperature. Units: K
            %
            temperature = calllib(ct, 'thermo_temperature', tp.tpID);
        end

        function pressure = get.P(tp)
            % GET.P  Get the pressure.
            % pressure = tp.P
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Pressure. Units: Pa
            %
            pressure = calllib(ct, 'thermo_pressure', tp.tpID);
        end

        function density = get.D(tp)
            % GET.D  Get the density.
            % density = tp.D
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Mass density. Units: kg/m**3
            %
            density = calllib(ct, 'thermo_density', tp.tpID);
        end

        function volume = get.V(tp)
            % GET.V  Get the specific volume.
            % volume = tp.V
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %    Specific volume depending on the basis. Units:
            %    m^3/kmol (molar) m^3/kg (mass).
            volume = 1/tp.D;
        end

        function moleFractions = get.X(tp)
            % GET.X  Get the mole fractions of all species.
            % moleFractions = tp.X
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Vector of species mole fractions for input phase. If
            %     no output argument is specified, a bar plot is produced.
            %
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'thermo_getMoleFractions', ...
                    tp.tpID, nsp, pt);
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
            % MOLEFRACTION  Get the mole fraction of one or a list of species.
            % x = tp.moleFraction(species)
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param species:
            %     String or cell array of strings of species whose mole
            %     fraction is desired
            % :return:
            %     Scalar or vector double mole fractions
            %
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
            % GET.Y  Get the mass fractions of all species.
            % massFractions = tp.Y
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Vector of species mass fractions for input phase. If
            %     no output argument is specified, a bar plot is produced.
            %
            nsp = tp.nSpecies;
            yy = zeros(1, nsp);
            pt = libpointer('doublePtr', yy);
            calllib(ct, 'thermo_getMassFractions', ...
                    tp.tpID, nsp, pt);
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
            % MASSFRACTION  Get the mass fraction of one or a list of species.
            % y = tp.massFraction(species)
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param species:
            %     String or cell array of strings of species whose mass
            %     fraction is desired
            % :return:
            %     Scalar or vector double mass fractions
            %
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

        function mu = chemicalPotentials(tp)
            % CHEMPOTENTIALS  Get the chemical potentials of the species.
            % mu = tp.chemPotentials
            % The expressions used to compute the chemical potential
            % depend on the model implemented by the underlying kernel
            % thermo manager.
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase).
            % :return:
            %     Vector of species chemical potentials. Units: J/kmol
            %
            %        This method returns an array containing the species
            %        chemical potentials [J/kmol]. The expressions used to
            %        compute these depend on the model implemented by the
            %        underlying kernel thermo manager.
            %
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'thermo_chemPotentials', ...
                    tp.tpID, nsp, pt);
            mu = pt.Value;
        end

        function c = cv(tp)
            % CV  Get the basis-dependent specific heat at constant volume.
            % c = tp.cv
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     specific heat of the mixture at constant volume.
            %     Units: J/kg-K (mass basis) or J/kmol-K (molar basis).
            %
            if strcmp(tp.basis, 'molar')
                c = calllib(ct, 'thermo_cv_mole', tp.tpID);
            elseif strcmp(tp.basis, 'mass')
                c = calllib(ct, 'thermo_cv_mass', tp.tpID);
            else error("basis not specified");
            end
        end

        function c = cp(tp)
            % CP  Get the basis-dependent specific heat at constant pressure.
            % v = tp.cp
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Specific heat of the mixture at constant pressure.
            %     Units: J/kg-K (mass-basis) or J/mol-K (molar-basis).
            %
            if strcmp(tp.basis, 'molar')
                c = calllib(ct, 'thermo_cp_mole', tp.tpID);
            elseif strcmp(tp.basis, 'mass')
                c = calllib(ct, 'thermo_cp_mass', tp.tpID);
            else error("basis not specified");
            end
        end

        function d = critDensity(tp)
            % CRITDENSITY  Get the critical density.
            % v = tp.critDensity
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Critical density. Units: kg/m**3
            %
            d = calllib(ct, 'thermo_critDensity', tp.tpID);
        end

        function p = critPressure(tp)
            % CRITPRESSURE  Get the critical pressure.
            % v = tp.critPressure
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Critical pressure. Units: Pa
            %
            p = calllib(ct, 'thermo_critPressure', tp.tpID);
        end

        function t = critTemperature(tp)
            % CRITTEMPERATURE  Get the critical temperature.
            % v = tp.critTemperature
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Critical temperature. Units: K
            %
            t = calllib(ct, 'thermo_critTemperature', tp.tpID);
        end

        function v = electricPotential(tp)
            % ELECTRICPOTENTIAL  Get the electric potential.
            % v = tp.electricPotential
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     The electric potential of the phase. Units: V
            %
            v = calllib(ct, 'thermo_electricPotential', tp.tpID);
        end

        function e = eosType(tp)
            % EOSTYPE  Get the type of the equation of state.
            % e = tp.eosType
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     An string identifying the equation of state.
            %
            buflen = calllib(ct, 'thermo_getEosType', tp.tpID, 0, '');
            if buflen > 0
                aa = char(zeros(1, buflen));
                [~, aa] = calllib(ct, 'thermo_getEosType', ...
                                        tp.tpID, buflen, aa);
            end
            e = aa;
        end

        function v = isIdealGas(tp)
            % ISIDEALGAS  Get a flag indicating whether the phase is an ideal gas.
            % v = tp.isIdealGas
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     True (1) if the phase is an ideal gas or ideal gas
            %     mixture, and false (0) otherwise.
            %
            if strcmp(tp.eosType, 'IdealGas')
                v = 1;
            else
                v = 0;
            end
            v = 1;
        end

        function b = isothermalCompressibility(tp)
            % ISOTHERMALCOMPRESSIBILITY  Get the isothermal compressibility.
            % b = tp.isothermalCompressibility
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Isothermal Compressibility. Units: 1/Pa
            %
            b = calllib(ct, 'thermo_isothermalCompressibility', tp.tpID);
        end

        function t = maxTemp(tp)
            % MAXTEMP  Get the maximum temperature of the parameter fits.
            % v = tp.maxTemp
            % The parameterizations used to represent the temperature-dependent
            % species thermodynamic properties are generally only valid in some
            % finite temperature range, which may be different for each species
            % in the phase.
            %
            % See also: :mat:func:`minTemp`
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Vector of maximum temperatures of all species
            %
            t = calllib(ct, 'thermo_maxTemp', tp.tpID, -1);
        end

        function t = minTemp(tp)
            % MINTEMP  Get the minimum temperature of the parameter fits.
            % v = tp.minTemp
            % The parameterizations used to represent the temperature-dependent
            % species thermodynamic properties are generally only valid in some
            % finite temperature range, which may be different for each species
            % in the phase.
            %
            % See also: :mat:func:`maxTemp`
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Vector of minimum temperatures of all species
            %
            t = calllib(ct, 'thermo_minTemp', tp.tpID, -1);
        end

        function p = refPressure(tp)
            % REFPRESSURE  Get the reference pressure.
            % v = tp.refPressure
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Reference pressure in Pa. All standard-state
            %     thermodynamic properties are for this pressure.
            %
            p = calllib(ct, 'thermo_refPressure', tp.tpID, -1);
        end

        function p = satPressure(tp, t)
            % SATPRESSURE  Get the saturation pressure for a given temperature.
            % v = tp.satPressure(T)
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param T:
            %     Temperature Units: K
            % :return:
            %     Saturation pressure for temperature T. Units: Pa
            %
            p = calllib(ct, 'thermo_satPressure', tp.tpID, t);
        end

        function t = satTemperature(tp, p)
            % SATTEMPERATURE  Get the saturation temperature for a given pressure.
            % v = tp.satTemperature(p)
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :param p:
            %     Pressure. Units: Pa
            % :return:
            %     Saturation temperature for pressure p. Units: K
            %
            t = calllib(ct, 'thermo_satTemperature', tp.tpID, p);
        end

        function c = soundspeed(tp)
            % SOUNDSPEED  Get the speed of sound.
            % c = tp.soundspeed
            % If the phase is an ideal gas, the speed of sound
            % is calculated by:
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
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :return:
            %     The speed of sound. Units: m/s
            %
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
            % THERMALEXPANSIONCOEFF  Get the thermal expansion coefficient.
            % a = tp.thermalExpansionCoeff
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :return:
            %     Thermal Expansion Coefficient. Units: 1/K
            %
            a = calllib(ct, 'thermo_thermalExpansionCoeff', tp.tpID);
        end

        function v = vaporFraction(tp)
            % VAPORFRACTION  Get the vapor fraction.
            % v = tp.vaporFraction
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :return:
            %     Vapor fraction.
            %
            v = calllib(ct, 'thermo_vaporFraction', tp.tpID);
        end

        function enthalpy = get.H(tp)
            % GET.H  Get the mass specific enthalpy.
            % enthalpy = tp.H
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %    Enthalpy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            if strcmp(tp.basis, 'molar')
                enthalpy = calllib(ct, 'thermo_enthalpy_mole', tp.tpID);
            elseif strcmp(tp.basis, 'mass')
                enthalpy = calllib(ct, 'thermo_enthalpy_mass', tp.tpID);
            else error("basis not specified");
            end
        end

        function enthalpy = enthalpies_RT(tp)
            % ENTHALPIES_RT  Get the non-dimensional enthalpies.
            % v = tp.enthalpies_RT
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Vector of standard-state species enthalpies
            %     divided by RT, where R is the universal gas
            %     constant and T is the temperature. For gaseous species, these
            %     values are ideal gas enthalpies.
            %
            nsp = tp.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'thermo_getEnthalpies_RT', tp.tpID, nsp, pt);
            enthalpy = pt.Value;
        end

        function entropy = get.S(tp)
            % GET.S  Get the mass specific entropy.
            % entropy = tp.S
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %     Entropy of the mixture depending on the basis.
            %     Units: J/kmol-K (molar) J/kg-K (mass)
            %
            if strcmp(tp.basis, 'molar')
                entropy = calllib(ct, 'thermo_entropy_mole', tp.tpID);
            elseif strcmp(tp.basis, 'mass')
                entropy = calllib(ct, 'thermo_entropy_mass', tp.tpID);
            else error("basis not specified");
            end
        end

        function intEnergy = get.U(tp)
            % GET.U  Get the mass specific internal energy.
            % intEnergy = tp.U
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %    Internal energy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            if strcmp(tp.basis, 'molar')
                intEnergy = calllib(ct, 'thermo_intEnergy_mole', tp.tpID);
            elseif strcmp(tp.basis, 'mass')
                intEnergy = calllib(ct, 'thermo_intEnergy_mass', tp.tpID);
            else error("basis not specified");
            end
        end

        function gibbs = get.G(tp)
            % GET.G  Get the mass specific Gibbs function.
            % gibbs = tp.G
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     object that derives from ThermoPhase)
            % :return:
            %    Gibbs free energy of the mixture depending on the basis.
            %    Units: J/kmol (molar) J/kg (mass).

            if strcmp(tp.basis, 'molar')
                gibbs = calllib(ct, 'thermo_gibbs_mole', tp.tpID);
            elseif strcmp(tp.basis, 'mass')
                gibbs = calllib(ct, 'thermo_gibbs_mass', tp.tpID);
            else error("basis not specified");
            end
        end

        %% ThermoGet multi methods

        function output = get.DP(tp)
            output = {tp.D, tp.P};
        end

        function output = get.DPX(tp)
            output = {tp.D, tp.P, tp.X};
        end

        function output = get.DPY(tp)
            output = {tp.D, tp.P, tp.Y};
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
            % SETELECTRICPOTENTIAL  Set the electric potential.
            % tp.setElectricPotential(phi)
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param phi:
            %     Electric potential. Units: V
            %
            calllib(ct, 'thermo_setElectricPotential', tp.tpID, phi);
        end

        function tp = setState_Psat(tp, p, q)
            % SETSTATEPSAT  Set the fluid state using the given pressure and quality.
            % tp.setState_Psat(p, q)
            % The fluid state will be set to a saturated liquid-vapor state using the
            % input pressure and vapor fraction (quality) as the independent,
            % intensive variables.
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param p:
            %     Pressure (Pa)
            % :param q:
            %     Vapor fraction
            %
            calllib(ct, 'thermo_setState_Psat', tp.tpID, p, q);
        end

        function tp = setState_Tsat(tp, t, q)
            % SETSTATETSAT  Set the fluid state using the given temperature and quality.
            % tp.setState_Tsat(t, q)
            % The fluid state will be set to a saturated liquid-vapor state using the
            % input temperature and vapor fraction (quality) as the independent,
            % intensive variables.
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param t:
            %     Temperature (K)
            % :param q:
            %     Vapor fraction
            %
            calllib(ct, 'thermo_setState_Tsat', tp.tpID, t, 1 - q);
        end

        function set.T(tp, temperature)
            % SET.T  Set the temperature.
            % tp.T = temperature
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param temperature:
            %     Temperature. Units: K
            %
            if temperature <= 0
                error('The temperature must be positive');
            end
            calllib(ct, 'thermo_setTemperature', tp.tpID, temperature);
        end

        function set.P(tp, pressure)
            % SET.P  Set the pressure.
            % tp.P = pressure
            % The pressure is set by changing the density holding the
            % temperature and chemical composition fixed.
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param pressure:
            %     Pressure. Units: Pa
            %
            if pressure <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_setPressure', tp.tpID, pressure);
        end

        function set.D(tp, density)
            % SET.D  Set the density.
            % tp.D = density
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param density:
            %     Density. Units: kg/m**3
            %
            if density <= 0
                error('The density must be positive');
            end
            calllib(ct, 'thermo_setDensity', tp.tpID, density);
        end

        function set.X(tp, xx)
            % SET.X  Set the species mole fractions.
            % tp.X = xx
            % Note that calling :mat:func:`setMoleFractions` leaves the temperature and
            % density unchanged, and therefore the pressure changes if the new
            % composition has a molar mass that is different than the old
            % composition. If it is desired to change the composition and hold
            % the pressure fixed, use method :mat:func:`set` and specify the mole
            % fractions and the pressure, or call :mat:func:`setPressure`
            % after calling :mat:func:`setmoleFractions`.
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param xx:
            %     Vector of mole fractions whose length must be the same as
            %     the number of species. Alternatively, a string in the format
            %     ``'SPEC:Y,SPEC2:Y2'`` that specifies the mole fraction of
            %     specific species.
            %
            tol = 1e-9;
            if isa(xx, 'double')
                nsp = tp.nSpecies;
                if sum(xx) - 1 <= tol
                    norm = 0;
                else norm = 1;
                end
                calllib(ct, 'thermo_setMoleFractions', tp.tpID, ...
                        nsp, xx, norm);
            elseif isa(xx, 'char')
                calllib(ct, 'thermo_setMoleFractionsByName', tp.tpID, xx);
            end
        end

        function set.Y(tp, yy)
            % SET.Y  Set the species mass fractions.
            % tp.Y = yy
            % Note that calling :mat:func:`setMassFractions` leaves the temperature and
            % density unchanged, and therefore the pressure changes if the new
            % composition has a molar mass that is different than the old
            % composition. If it is desired to change the composition and hold
            % the pressure fixed, use method :mat:func:`set` and specify the mass
            % fractions and the pressure, or call :mat:func:`setPressure`
            % after calling :mat:func:`setMassFractions`.
            %
            % :param tp:
            %     Instance of class :mat:func:`ThermoPhase` (or another
            %     class derived from ThermoPhase)
            % :param yy:
            %     Vector of mass fractions whose length must be the same as
            %     the number of species. Alternatively, a string in the format
            %     ``'SPEC:Y,SPEC2:Y2'`` that specifies the mass fraction of
            %     specific species.
            %
            if isa(yy, 'double')
                nsp = tp.nSpecies;
                if sum(yy) -1 <= 1e-9
                    norm = 0;
                else norm = 1;
                end
                calllib(ct, 'thermo_setMassFractions', tp.tpID, ...
                        nsp, yy, norm);
            elseif isa(yy, 'char')
                calllib(ct, 'thermo_setMassFractionsByName', tp.tpID, yy);
            end
        end

        %% PhaseSet multi methods

        function set.DP(tp, input)
            d = input{1};
            p = input{2};
            if d <= 0
                error('The density must be positive');
            end
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_RP', tp.tpID, [d, p]);
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
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_HP', tp.tpID, [h, p]);
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
            if p <= 0
                error('The pressure must be positive');
            end
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_PV', tp.tpID, [p, v]);
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
            s = input{1};
            h = input{2};
            calllib(ct, 'thermo_set_SH', tp.tpID, [s, h]);
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
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_SP', tp.tpID, [s, p]);
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
            if t <= 0
                error('The temperature must be positive');
            end
            calllib(ct, 'thermo_set_ST', tp.tpID, [s, t]);
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
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_SV', tp.tpID, [s, v]);
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
            t = input{1};
            if t <= 0
                error('The temperature must be positive');
            end
            h = input{2};
            calllib(ct, 'thermo_set_TH', tp.tpID, [t, h]);
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
            t = input{1};
            v = input{2};
            if t <= 0
                error('The temperature must be positive');
            end
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_TV', tp.tpID, [t, v]);
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
            if p <= 0
                error('The pressure must be positive');
            end
            calllib(ct, 'thermo_set_UP', tp.tpID, [u, p]);
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
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_UV', tp.tpID, [u, v]);
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
            if v <= 0
                error('The specific volume must be positive');
            end
            calllib(ct, 'thermo_set_VH', tp.tpID, [v, h]);
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

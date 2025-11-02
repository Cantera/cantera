classdef (Abstract) ReactorBase < handle

    properties (SetAccess = immutable)

        id % ID of reactor.

        % The Solution object used to represent the contents of this reactor.
        % If the `clone` parameter of the reactor constructor is true,
        % a new instance of :mat:class:`Solution` will be created for this reactor.
        % Otherwise, the same instance as passed to the constructor will be used.
        phase

    end

    properties (SetAccess = public)

        type  % reactor type.

        name  % Name of reactor.

        % Density of the reactor contents at the end of the last call to
        % "advance" or "step". Unit: kg/m^3.
        D

        % Pressure of the reactor contents at the end of the last call to
        % "advance" or "step". Unit: K.
        P

        % Mass of the reactor contents at the end of the last call to
        % "advance" or "step". Unit: kg.
        M

        % Mass specific enthalpy of the reactor contents at the end of the last call to
        % "advance" or "step". Unit: J/kg.
        H

        % Mass specific internal energy of the reactor contents at the end of the last call to "advance" or "step". Unit: J/kg.
        U

        % Temperature of the reactor contents at the end of the last call to
        % "advance" or "step". Unit: K.
        T

        % Volume of the reactor contents at the end of the last call to
        % "advance" or "step". Unit: m^3.
        V

        % Mass fractions of the reactor contents at the end of the last call to
        % "advance" or "step".
        Y

        % Enable or disable changing reactor composition by reactions. ::
        %
        %     >> r.chemistryEnabled = flag
        %
        % If the chemistry is disabled, no reactions take place. Note that the
        % composition can still change due to mass transfer at inlets and outlets.
        %
        % :param r:
        %    Instance of :mat:class:`ReactorBase`.
        % :param flag:
        %    Boolean to enable or disable chemical reactions.
        chemistryEnabled

        % Enable or disable solving the energy equation. ::
        %
        %     >> r.energyEnabled = flag
        %
        % If the energy equation is disabled, then the reactor temperature is constant.
        %
        % :param r:
        %    Instance of :mat:class:`ReactorBase`.
        % :param flag:
        %    Boolean to enable or disable energy equations.
        energyEnabled

        area % Area of the reactor in m².

    end

    methods
        %% ReactorBase Class Constructor

        function r = ReactorBase(id)
            % ReactorBase Class ::
            %
            % :param id:
            %    ID of reactor object instantiated by specialization.

            arguments
                id (1,1) {mustBeInteger, mustBeGreaterThan(id, -1)}
            end

            r.id = id;
            phaseID = ctFunc('mReactor_phase', id);
            r.phase= Solution(phaseID);
        end

        %% ReactorBase Class Destructor

        function delete(r)
            % Delete the :mat:class:`ReactorBase` object.
            ctFunc('mReactor_del', r.id);
        end

        %% ReactorBase Utility Methods

        function addSensitivityReaction(r, m)
            % Specifies that the sensitivity of the state variables with
            % respect to reaction m should be computed. The reactor must be
            % part of a network first. ::
            %
            %     >> r.addSensitivityReaction(m)
            %
            % :param m:
            %    Index number of reaction.

            ctFunc('mReactor_addSensitivityReaction', r.id, m);
        end

        %% ReactorBase Get Methods

        function typ = get.type(r)
            typ = ctString('mReactor_type', r.id);
        end

        function name = get.name(r)
            name = ctString('mReactor_name', r.id);
        end

        function temperature = get.T(r)
            temperature = ctFunc('mReactor_temperature', r.id);
        end

        function pressure = get.P(r)
            pressure = ctFunc('mReactor_pressure', r.id);
        end

        function rho = get.D(r)
            rho = ctFunc('mReactor_density', r.id);
        end

        function mass = get.M(r)
            mass = ctFunc('mReactor_mass', r.id);
        end

        function volume = get.V(r)
            volume = ctFunc('mReactor_volume', r.id);
        end

        function enthalpy_mass = get.H(r)
            enthalpy_mass = ctFunc('mReactor_enthalpy_mass', r.id);
        end

        function intEnergy_mass = get.U(r)
            intEnergy_mass = ctFunc('mReactor_intEnergy_mass', r.id);
        end

        function a = get.area(r)
            a = ctFunc('mReactor_area', r.id);
        end

        function yi = massFraction(r, species)
            % Get the mass fraction of a species.
            %
            % :param species:
            %    String/Char array name of the species.

            k = r.phase.speciesIndex(species) - 1;
            yi = ctFunc('mReactor_massFraction', r.id, k);
        end

        function massFractions = get.Y(r)
            massFractions = ctArray('mReactor_massFractions', r.phase.nSpecies, r.id);
        end

        function flag = get.chemistryEnabled(r)
            flag = logical(ctFunc('mReactor_chemistryEnabled', r.id));
        end

        function flag = get.energyEnabled(r)
            flag = logical(ctFunc('mReactor_energyEnabled', r.id));
        end

        %% ReactorBase set methods

        function set.name(r, name)
            ctFunc('mReactor_setName', r.id, name);
        end

        function set.V(r, v0)
            ctFunc('mReactor_setInitialVolume', r.id, v0);
        end

        function set.chemistryEnabled(r, flag)
            arguments
                r {mustBeA(r, 'ReactorBase')}
                flag (1,1) logical
            end

            ctFunc('mReactor_setChemistryEnabled', r.id, flag);
        end

        function set.energyEnabled(r, flag)
            arguments
                r {mustBeA(r, 'ReactorBase')}
                flag (1,1) logical
            end

            ctFunc('mReactor_setEnergyEnabled', r.id, flag);
        end

        function set.area(s, a)
            ctFunc('mReactor_setArea', s.id, a);
        end

    end

end

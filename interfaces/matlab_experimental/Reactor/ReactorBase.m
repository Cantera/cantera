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
        %     >> r.chemistry = flag
        %
        % If the chemistry is disabled, then the reactor composition is
        % constant. The parameter should be the string "on" to enable
        % the species equations, or "off" to disable it.
        %
        % By default, :mat:class:`ReactorBase` objects are created with the species
        % equations enabled if there are reactions present in the
        % mechanism file, and disabled otherwise.
        %
        % :param r:
        %    Instance of :mat:class:`ReactorBase`.
        % :param flag:
        %    String, either "on" or "off" to enable or disable chemical
        %    reactions, respectively.
        chemistry

        % Enable or disable solving the energy equation. ::
        %
        %     >> r.energy = flag
        %
        % If the energy equation is disabled, then the reactor
        % temperature is constant. The parameter should be the string
        % "on" to enable the energy equation, or "off" to disable it.
        %
        % By default, :mat:class:`ReactorBase` objects are created with the energy
        % equation enabled, so usually this method.
        %
        % :param r:
        %    Instance of :mat:class:`ReactorBase`.
        % :param flag:
        %    String, either "on" or "off" to enable or disable chemical
        %    reactions, respectively.
        energy

        massFlowRate % Mass flow rate in kg/s.

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
            phaseID = ctFunc('reactor_phase', id);
            r.phase= Solution(phaseID);
        end

        %% ReactorBase Class Destructor

        function delete(r)
            % Delete the :mat:class:`ReactorBase` object.
            ctFunc('reactor_del', r.id);
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

            ctFunc('reactor_addSensitivityReaction', r.id, m);
        end

        %% ReactorBase Get Methods

        function typ = get.type(r)
            typ = ctString('reactor_type', r.id);
        end

        function name = get.name(r)
            name = ctString('reactor_name', r.id);
        end

        function temperature = get.T(r)
            temperature = ctFunc('reactor_temperature', r.id);
        end

        function pressure = get.P(r)
            pressure = ctFunc('reactor_pressure', r.id);
        end

        function rho = get.D(r)
            rho = ctFunc('reactor_density', r.id);
        end

        function mass = get.M(r)
            mass = ctFunc('reactor_mass', r.id);
        end

        function volume = get.V(r)
            volume = ctFunc('reactor_volume', r.id);
        end

        function enthalpy_mass = get.H(r)
            enthalpy_mass = ctFunc('reactor_enthalpy_mass', r.id);
        end

        function intEnergy_mass = get.U(r)
            intEnergy_mass = ctFunc('reactor_intEnergy_mass', r.id);
        end

        function yi = massFraction(r, species)
            % Get the mass fraction of a species.
            %
            % :param species:
            %    String or one-based integer id of the species.

            if ischar(species)
                k = r.phase.speciesIndex(species) - 1;
            else
                k = species - 1;
            end

            yi = ctFunc('reactor_massFraction', r.id, k);
        end

        function massFractions = get.Y(r)

            nsp = r.phase.nSpecies;
            massFractions = zeros(1, nsp);

            for i = 1:nsp
                massFractions(i) = r.massFraction(i);
            end

        end

        %% ReactorBase set methods

        function set.name(r, name)
            ctFunc('reactor_setName', r.id, name);
        end

        function set.V(r, v0)
            ctFunc('reactor_setInitialVolume', r.id, v0);
        end

        function set.massFlowRate(r, MFR)
            ctFunc('reactor_setMassFlowRate', r.id, MFR);
            r.massFlowRate = MFR;
        end

        function set.chemistry(r, flag)

            if strcmp(flag, 'on')
                cflag = true;
            elseif strcmp(flag, 'off')
                cflag = false;
            else
                error('Input must be "on" or "off"');
            end

            ctFunc('reactor_setChemistry', r.id, cflag);
        end

        function set.energy(r, flag)

            if strcmp(flag, 'on')
                eflag = 1;
            elseif strcmp(flag, 'off')
                eflag = 0;
            else
                error('Input must be "on" or "off".');
            end

            ctFunc('reactor_setEnergy', r.id, eflag);

        end

    end

end

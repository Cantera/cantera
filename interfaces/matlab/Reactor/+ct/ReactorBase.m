classdef (Abstract) ReactorBase < handle

    properties (SetAccess = immutable)

        id = -1  % ID of reactor.

        % The Solution object used to represent the contents of this reactor.
        % If the `clone` parameter of the reactor constructor is true,
        % a new instance of :mat:class:`Solution` will be created for this reactor.
        % Otherwise, the same instance as passed to the constructor will be used.
        phase

    end

    properties (SetAccess = public)

        type  % reactor type.

        name  % Name of reactor.

        % Density [kg/m³] of the reactor contents at the end of the last call to
        % "advance" or "step".
        D

        % Pressure [Pa] of the reactor contents at the end of the last call to
        % "advance" or "step".
        P

        % Mass [kg] of the reactor contents at the end of the last call to
        % "advance" or "step".
        M

        % Mass specific enthalpy [J/kg] of the reactor contents at the end of the last
        % call to "advance" or "step".
        H

        % Mass specific internal energy [J/kg] of the reactor contents at the end of the
        % last call to "advance" or "step".
        U

        % Temperature [K] of the reactor contents at the end of the last call to
        % "advance" or "step".
        T

        % Volume [m³] of the reactor contents at the end of the last call to "advance"
        % or "step".
        V

        % Mass fractions of the reactor contents at the end of the last call to
        % "advance" or "step".
        Y

        % Boolean to enable or disable changing reactor composition by reactions.
        %
        % If the chemistry is disabled, no reactions take place. Note that the
        % composition can still change due to mass transfer at inlets and outlets.
        chemistryEnabled

        % Boolean to enable or disable solving the energy equation.
        %
        % If the energy equation is disabled, then the reactor temperature is constant.
        energyEnabled

        area % Area of the reactor [m²].

    end

    methods
        %% ReactorBase Class Constructor

        function r = ReactorBase(id)
            % ReactorBase Class
            %
            % :param id:
            %    ID of reactor object instantiated by specialization.

            arguments
                id (1,1) {mustBeInteger, mustBeGreaterThan(id, -1)}
            end

            r.id = id;
            phaseID = ct.impl.call('mReactor_phase', id);
            r.phase = ct.Solution(phaseID);
        end

        %% ReactorBase Class Destructor

        function delete(obj)
            % Delete the :mat:class:`ReactorBase` object.
            if obj.id >= 0
                ct.impl.call('mReactor_del', obj.id);
            end
        end

        %% ReactorBase Utility Methods

        function addSensitivityReaction(obj, m)
            % Specifies that the sensitivity of the state variables with
            % respect to reaction m should be computed. The reactor must be
            % part of a network first. ::
            %
            %     >> r.addSensitivityReaction(m)
            %
            % :param m:
            %    Index number of reaction.

            ct.impl.call('mReactor_addSensitivityReaction', obj.id, m);
        end

        %% ReactorBase Get Methods

        function typ = get.type(obj)
            typ = ct.impl.getString('mReactor_type', obj.id);
        end

        function name = get.name(obj)
            name = ct.impl.getString('mReactor_name', obj.id);
        end

        function temperature = get.T(obj)
            temperature = ct.impl.call('mReactor_temperature', obj.id);
        end

        function pressure = get.P(obj)
            pressure = ct.impl.call('mReactor_pressure', obj.id);
        end

        function rho = get.D(obj)
            rho = ct.impl.call('mReactor_density', obj.id);
        end

        function mass = get.M(obj)
            mass = ct.impl.call('mReactor_mass', obj.id);
        end

        function volume = get.V(obj)
            volume = ct.impl.call('mReactor_volume', obj.id);
        end

        function enthalpy_mass = get.H(obj)
            enthalpy_mass = ct.impl.call('mReactor_enthalpy_mass', obj.id);
        end

        function intEnergy_mass = get.U(obj)
            intEnergy_mass = ct.impl.call('mReactor_intEnergy_mass', obj.id);
        end

        function a = get.area(obj)
            a = ct.impl.call('mReactor_area', obj.id);
        end

        function yi = massFraction(obj, species)
            % Get the mass fraction of a species.
            %
            % :param species:
            %    String/Char array name of the species.

            k = obj.phase.speciesIndex(species) - 1;
            yi = ct.impl.call('mReactor_massFraction', obj.id, k);
        end

        function massFractions = get.Y(obj)
            massFractions = ct.impl.getArray('mReactor_massFractions', ...
                                             obj.phase.nSpecies, obj.id);
        end

        function flag = get.chemistryEnabled(obj)
            flag = logical(ct.impl.call('mReactor_chemistryEnabled', obj.id));
        end

        function flag = get.energyEnabled(obj)
            flag = logical(ct.impl.call('mReactor_energyEnabled', obj.id));
        end

        %% ReactorBase set methods

        function set.name(obj, name)
            ct.impl.call('mReactor_setName', obj.id, name);
        end

        function set.V(obj, v0)
            ct.impl.call('mReactor_setInitialVolume', obj.id, v0);
        end

        function set.chemistryEnabled(obj, flag)
            arguments
                obj (1,1) ct.ReactorBase
                flag (1,1) logical
            end

            ct.impl.call('mReactor_setChemistryEnabled', obj.id, flag);
        end

        function set.energyEnabled(obj, flag)
            arguments
                obj (1,1) ct.ReactorBase
                flag (1,1) logical
            end

            ct.impl.call('mReactor_setEnergyEnabled', obj.id, flag);
        end

        function set.area(obj, a)
            ct.impl.call('mReactor_setArea', obj.id, a);
        end

    end

end

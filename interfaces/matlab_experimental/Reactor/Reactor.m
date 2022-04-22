classdef Reactor < handle

    properties
        type
        id
        contents
        D
        P
        M
        H
        U
        T
        V
        Y
        Mdot
        ChemistryFlag
        EnergyFlag
    end

    methods
        %% Reactor class constructor

        function r = Reactor(content, typ)
            % A 'Reactor' object simulates a perfectly-stirred reactor. It
            % has a time-dependent tstate, and may be coupled to other
            % reactors through flow lines or through walls that may expand
            % or contract and/orconduct heat.
            %
            % :parameter contents:
            %    Instance of class 'Solution' representing the contents of
            %    the reactor.
            % :parameter typ:
            %    Character array of reactor type. Options are:
            %    'Reservoir'
            %    'Reactor'
            %    'FlowReactor'
            %    'ConstPressureReactor'
            %    'IdealGasReactor'
            %    'IdealGasConstPressureReactor'
            % :return:
            %    Instance of class 'Reactor'.

            checklib;

            if nargin == 0
                content = 0;
                typ = 'Reactor';
            elseif nargin == 1
                typ = 'Reactor'
            elseif nargin > 2
                error('too many arguments');
            end

            r.type = char(typ);
            r.id = calllib(ct, 'reactor_new', typ);

%             if r.id < 0
%                 error(geterr);
%             end

            if isa(content, 'Solution')
                r.insert(content);
            elseif ~(isa(contents, 'double') && contents == 0)
                error('Reactor contents must be an object of type "Solution"');
            end

        end

        %% Utility methods

        function clear(r)
            % Clear the reactor from memory.

            calllib(ct, 'reactor_del', r.id);
        end

        function addSensitivityReaction(r, m)
            % Specifies that the sensitivity of the state variables with
            % respect to reaction m should be computed. The reactor must be
            % part of a network first.
            %
            % :parameter m:
            %    Index number of reaction.

            calllib(ct, 'reactor_addSensitivityReaction', r.id, m);
        end

        function insert(r, gas)
            % Insert a solution or mixture into a reactor.
            %
            % :parameter r:
            %    Instance of class 'Reactor'.
            % :parameter gas:
            %    Instance of class 'Solution'.

            r.contents = gas;
            setThermoMgr(r, gas);
            if ~strcmp(r.type, 'Reservoir')
                setKineticsMgr(r, gas);
            end
        end

        function setThermoMgr(r, t)
            % Set the thermodynamics manager.
            %
            % This method is used internally during Reactor initialization,
            % but is usually not called by users.
            %
            % :parameter r:
            %    Instance of class 'Reactor'.
            % :parameter t:
            %    Instance of class 'ThermoPhase' or another object
            %    containing an instance of that class.

            if ~isa(t, 'ThermoPhase')
                error('Wrong object type');
            end

            calllib(ct, 'reactor_setThermoMgr', r.id, t.tpID);
        end

        function setKineticsMgr(r, k)
            % Set the kinetics manager.
            %
            % This method is used internally during Reactor initialization,
            % but is usually not called by users.
            %
            % :parameter r:
            %    Instance of class 'Reactor'.
            % :parameter t:
            %    Instance of class 'Kinetics' or another object
            %    containing an instance of that class.

            if ~isa(k, 'Kinetics')
                error('Wrong object type');
            end

            calllib(ct, 'reactor_setKineticsMgr', r.id, k.kinID);
        end

        %% Reactor get methods

        function temperature = get.T(r)
            % Get the temperature of the reactor.
            %
            % :return:
            %    The temperature of the reactor contents at the end of the
            %    last call to 'advance' or 'step'. Unit: K.

            temperature = calllib(ct, 'reactor_temperature', r.id);
        end

        function pressure = get.P(r)
            % Get the pressure of the reactor.
            %
            % :return:
            %    The pressure of the reactor contents at the end of the
            %    last call to 'advance' or 'step'. Unit: Pa.

            pressure = calllib(ct, 'reactor_pressure', r.id);
        end

        function rho = get.D(r)
            % Get the density of the reactor.
            %
            % :return:
            %    Density of the phase in the input. Unit: kg/m^3.

            rho = calllib(ct, 'reactor_density', r.id);
        end

        function mass = get.M(r)
            % Get the mass of the reactor.
            %
            % :return:
            %    The mass of the reactor contents at the end of the
            %    last call to 'advance' or 'step'. The mass is retrieved
            %    from the solution vector. Unit: kg.

            mass = calllib(ct, 'reactor_mass', r.id);
        end

        function volume = get.V(r)
            % Get the volume of the reactor.
            %
            % :return:
            %    The volume of the reactor contents at the end of the
            %    last call to 'advance' or 'step'. Unit: m^3.

            volume = calllib(ct, 'reactor_volume', r.id);
        end

        function enthalpy_mass = get.H(r)
            % Get the mass specific enthalpy of the reactor.
            %
            % :return:
            %    The mass specific enthalpy of the reactor contents at the
            %    end of the last call to 'advance' or 'step'. The enthalpy
            %    is retrieved from the solution vector. Unit: J/kg.

            enthalpy_mass = calllib(ct, 'reactor_enthalpy_mass', r.id);
        end

        function intEnergy_mass = get.U(r)
            % Get the mass specific internal energy of the reactor.
            %
            % :return:
            %    The mass specific internal energy of the reactor contents
            %    at the end of the last call to 'advance' or 'step'. The
            %    internal energy is retrieved from the solution vector.
            %    Unit: J/kg.

            intEnergy_mass = calllib(ct, 'reactor_intEnergy_mass', r.id);
        end

        function yi = massFraction(r, species)
            % Get the mass fraction of a species.
            %
            % :parameter species:
            %    String or one-based integer id of the species.

            if ischar(species)
                k = r.contents.speciesIndex(species) - 1;
            else k = species - 1;
            end

            yi = calllib(ct, 'reactor_massFraction', r.id, k);
        end

        function massFractions = get.Y(r)
            % Get the mass fractions of the reactor.
            %
            % :return:
            %    The mass fractions of the reactor contents at the end of
            %    the last call to 'advance' or 'step'.

            nsp = r.contents.nSpecies;
            massFractions = zeros(1, nsp);
            for i = 1:nsp
                massFractions(i) = r.massFraction(i);
            end
        end

        %% Reactor set methods

        function setInitialVolume(r, v0)
            % Set the initial reactor volume.
            %
            % :parameter v0:
            %    Initial volume in m^3.

            calllib(ct, 'reactor_setInitialVolume', r.id, v0);
        end

        function r = set.Mdot(r, MFR)
            % Set the mass flow rate.
            %
            % :parameter MFR:
            %    Mass flow rate in kg/s.

            calllib(ct, 'reactor_setMassFlowRate', r.id, MFR);
            r.Mdot = MFR;
        end

        function r = set.ChemistryFlag(r, flag)
            % Enable or disable changing reactor composition by reactions.
            %
            % If the chemistry is disabled, then the reactor composition is
            % constant. The :parameter should be the string "on" to enable
            % the species equaionts, or "off" to disable it.
            %
            % By default, Reactor objects are created with the species
            % equations enabled if there are reactions present in the
            % mechanism file, and disabled otherwise.
            %
            % :parameter r:
            %    Instance of class 'Reactor'.
            % :parameter flag:
            %    String, either "on" or "off" to enable or disable chemical
            %    reactions, respectively.

            if strcmp(flag, 'on')
                iflag = true;
            elseif strcmp(flag, 'off')
                iflag = false;
            else error('Input must be "on" or "off"');
            end

            calllib(ct, 'reactor_setChemistry', r.id, iflag);
            r.ChemistryFlag = flag;
        end

        function r = set.EnergyFlag(r, flag)
            % Enable or disable solving the energy equation.
            %
            % If the energy equation is disabled, then the reactor
            % temperature is constant. The :parameter should be the string
            % "on" to enable the energy equation, or "off" to disable it.
            %
            % By default, Reactor objects are created with the energy
            % equation enabled, so usually this method
            %
            % :parameter r:
            %    Instance of class 'Reactor'.
            % :parameter flag:
            %    String, either "on" or "off" to enable or disable chemical
            %    reactions, respectively.

            iflag = -1;
            if strcmp(flag, 'on')
                iflag = 1;
            elseif strcmp(flag, 'off')
                iflag = 0;
            end

            if iflag >= 0
                calllib(ct, 'reactor_setEnergy', r.id, iflag);
            else error('Input must be "on" or "off".');
            end

            r.EnergyFlag = flag;
        end

    end
end

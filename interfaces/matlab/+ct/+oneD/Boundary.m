classdef (Abstract) Boundary < ct.oneD.Domain
    % Boundary Class.
    %
    % Base class for objects representing domain boundaries. The constructor is called
    % by derived classes and cannot be used directly.
    %
    % :param type:
    %    String type of Boundary. Possible values are:
    %      - `inlet`
    %      - `outlet`
    %      - `reacting-surface`
    %      - `surface`
    %      - `symmetry-plane`
    %      - `outlet-reservoir`
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution` or :mat:class:`ct.Interface`.
    % :param name:
    %     String, ID of the flow.

    properties
        massFlux % The mass flux [kg/s/m²] in the domain.

        T % Boundary temperature [K].

    end

    methods

        %% Boundary Class Constructor

        function obj = Boundary(type, phase, name)
            arguments
                type (1,1) string
                phase (1,1) ct.Solution
                name (1,1) string
            end
            id = ct.impl.call('mDomain_newBoundary1D', type, phase.solnID, name);
            obj@ct.oneD.Domain(id);

        end

        %% Boundary Methods

        function temperature = get.T(obj)
            temperature = ct.impl.call('mBdry_temperature', obj.domainID);
        end

        function set.T(obj, t)
            ct.impl.call('mBdry_setTemperature', obj.domainID, t);
        end

        function mdot = get.massFlux(obj)
            mdot = ct.impl.call('mBdry_mdot', obj.domainID);
        end

        function set.massFlux(obj, mdot)
            ct.impl.call('mBdry_setMdot', obj.domainID, mdot);
        end

        function y = massFraction(obj, k)
            % Get the mass fraction of a species given its integer index. ::
            %
            %     >> y = d.massFraction(k)
            %
            % This method returns the mass fraction of species ``k``, where
            % k is the integer index of the species in the flow domain
            % to which the boundary domain is attached.
            %
            % :param k:
            %     Integer species index.
            % :return:
            %     Mass fraction of species.

            y = ct.impl.call('mBdry_massFraction', obj.domainID, k - 1);

        end

        function v = value(obj, component)
            % Get the value of a component at a boundary. ::
            %
            %     >> d.value(component)
            %
            % :param component:
            %    String component for which the solution is desired.
            % :return:
            %    Value of the component in domain d.

            v = ct.impl.call('mDomain_value', obj.domainID, component);

        end

    end

end

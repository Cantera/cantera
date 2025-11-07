classdef (Abstract) Boundary1D < Domain1D
    % Create a Boundary domain. ::
    %
    %     >> m = Boundary(type, phase, name)
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
    %     Instance of class :mat:class:`Solution` or :mat:class:`Interface`.
    % :param name:
    %     String, ID of the flow.

    properties
        massFlux % The mass flux [kg/s/m²] in the domain.

        T % Boundary temperature [K].

    end

    methods

        %% Boundary Class Constructor

        function b = Boundary1D(type, phase, name)
            arguments
                type (1,1) string
                phase (1,1) Solution
                name (1,1) string
            end
            id = ctFunc('mDomain_newBoundary1D', type, phase.solnID, name);
            b@Domain1D(id);

        end

        %% Boundary Methods

        function temperature = get.T(obj)
            temperature = ctFunc('mBdry_temperature', obj.domainID);
        end

        function set.T(obj, t)
            ctFunc('mBdry_setTemperature', obj.domainID, t);
        end

        function mdot = get.massFlux(obj)
            mdot = ctFunc('mBdry_mdot', obj.domainID);
        end

        function set.massFlux(obj, mdot)
            ctFunc('mBdry_setMdot', obj.domainID, mdot);
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

            y = ctFunc('mBdry_massFraction', obj.domainID, k - 1);

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

            v = ctFunc('mDomain_value', obj.domainID, component);

        end

    end

end

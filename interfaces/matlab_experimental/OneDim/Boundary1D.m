classdef Boundary1D < Domain1D
    % Create a Boundary domain. ::
    %
    %     >> m = Boundary(type, phase, id)
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
    % :param id:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`Boundary`.

    properties

        % Mass flux. ::
        %
        %     >> mdot = d.massFlux
        %
        % :param d:
        %     Instance of class :mat:class:`Boundary`.
        % :return:
        %     The mass flux in the domain. Unit: [kg/s/m^2].
        massFlux

        T % Boundary Temperature. Units: K.

    end

    methods

        %% Boundary Class Constructor

        function b = Boundary1D(type, phase, id)

            ctIsLoaded;
            domainID = ctFunc('domain_newBoundary1D', type, phase.solnID, id);
            b@Domain1D(domainID);

        end

        %% Boundary Methods

        function temperature = get.T(d)
            temperature = ctFunc('bdry_temperature', d.domainID);
        end

        function set.T(d, t)
            ctFunc('bdry_setTemperature', d.domainID, t);
        end

        function mdot = get.massFlux(d)
            mdot = ctFunc('bdry_mdot', d.domainID);
        end

        function set.massFlux(d, mdot)
            ctFunc('bdry_setMdot', d.domainID, mdot);
        end

        function y = massFraction(d, k)
            % Get the mass fraction of a species given its integer index. ::
            %
            %     >> y = d.massFraction(k)
            %
            % This method returns the mass fraction of species ``k``, where
            % k is the integer index of the species in the flow domain
            % to which the boundary domain is attached.
            %
            % :param d:
            %     Instance of class :mat:class:`Boundary`.
            % :param k:
            %     Integer species index.
            % :return:
            %     Mass fraction of species.

            if d.domainIndex == 0
                error('No flow domain attached')
            end

            if strcmp(d.domainType, 'inlet')
                y = ctFunc('bdry_massFraction', d.domainID, k - 1);
            else
                error('Input domain must be an inlet');
            end

        end

        function setMoleFractions(d, x)
            % Set the mole fractions. ::
            %
            %     >> d.setMoleFractions(x)
            %
            % :param d:
            %     Instance of class :mat:class:`Boundary`.
            % :param x:
            %     String specifying the species and mole fractions in
            %     the format ``'SPEC:X,SPEC2:X2'``.

            ctFunc('bdry_setMoleFractions', d.domainID, x);
        end

    end

end

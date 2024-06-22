classdef Flow1D < Domain1D
    % Create a Flow domain. ::
    %
    %     >> m = Flow1D(type, phase, id)
    %
    % :param type:
    %    String type of domain. Possible values are:
    %      - `axisymmetric-flow`
    %      - `free-flow`
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`Flow1D`.

    properties
        P % Flow Pressure. Units: Pa.
    end

    methods

        %% Flow1D Class Constructor

        function f = Flow1D(type, phase, id)

            f@Domain1D(type, phase, id);

        end

        %% Flow1D Class Methods

        function pressure = get.P(d)
            pressure = ctFunc('flow1D_pressure', d.domainID);
        end

        function set.P(d, p)
            ctFunc('flow1D_setPressure', d.domainID, p);
        end

        function setFixedTempProfile(d, profile)
            % Set a fixed temperature profile. ::
            %
            %     >> d.setFixedTempProfile(profile)
            %
            % Set the temperature profile to use when the energy equation
            % is not being solved. The profile must be entered as an array of
            % positions / temperatures, which may be in rows or columns.
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param profile:
            %     :math:`n\times 2` or :math:`2\times n` array of ``n`` points
            %     at which the temperature is specified.

            sz = size(profile);

            if sz(1) == 2
                l = length(profile(1, :));
                ctFunc('flow1D_setFixedTempProfile', d.domainID, ...
                        l, profile(1, :), l, profile(2, :));
            elseif sz(2) == 2
                l = length(profile(:, 1));
                ctFunc('flow1D_setFixedTempProfile', d.domainID, ...
                        l, profile(:, 1), l, profile(:, 2));
            else
                error('Wrong temperature profile array shape.');
            end

        end

    end

end

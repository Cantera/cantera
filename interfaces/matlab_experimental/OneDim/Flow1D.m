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

            ctIsLoaded;
            domainID = ctFunc('domain_newFlow1D', type, phase.solnID, id);
            f@Domain1D(domainID);
            f.energyEnabled = false;
            f.soretEnabled = false;

        end

        %% Flow1D Class Methods

        function pressure = get.P(d)
            pressure = ctFunc('flow_pressure', d.domainID);
        end

        function set.P(d, p)
            ctFunc('flow_setPressure', d.domainID, p);
        end

        function v = values(d, component)
            % Get the value of a component in a domain. ::
            %
            %     >> d.values(component)
            %
            % :param d:
            %    Instance of class :mat:class:`Flow1D`.
            % :param component:
            %    String component for which the solution is desired.
            % :return:
            %    Value of the component in domain d.

            v = ctArray('domain_values', d.domainID, component);

        end

        function setProfile(d, component, zProfile, vProfile)
            % Specify a profile for one component. ::
            %
            %     >> d.setProfile(component, zProfile, vProfile)
            %
            % The solution vector values for this component will be linearly
            % interpolated from the discrete function defined by two input arrays.
            % This method can be called at any time, but is usually used to set the
            % initial guess for the solution.
            %
            % Example (assuming 'd' is an instance of class :mat:class:`Domain1D`):
            %    >> zr = [0.0, 0.1, 0.2, 0.4, 0.8, 1.0];
            %
            %    >> v = [500, 650, 700, 730, 800, 900];
            %
            %    >> d.setProfile('T', zr, v);
            %
            % :param d:
            %    Instance of class :mat:class:`Domain1D`.
            % :param component:
            %    Component name.
            % :param zProfile:
            %    Relative (normalized) positions. "zProfile(1) = 0.0" corresponds to the
            %    leftmost grid point in the specified domain, and "zProfile(n) = 1.0"
            %    corresponds to the rightmost grid point
            % :param vProfile:
            %    Array containing component values at positions.

            ctFunc('domain_setProfile', d.domainID, ...
                   component, zProfile, vProfile);

        end

        function setFlatProfile(d, component, v)
            % Set a component to a value across the entire domain. ::
            %
            %     >> d.setFlatProfile(component, v)
            %
            % :param d:
            %    Instance of class :mat:class:`Domain1D`.
            % :param component:
            %    Component to be set.
            % :param v:
            %    Double value to be set.

            ctFunc('domain_setFlatProfile', d.domainID, component, v);
        end

        function setFixedTempProfile(d, zFixed, tFixed)
            % Set a fixed temperature profile. ::
            %
            %     >> d.setFixedTempProfile(zFixed, tFixed)
            %
            % Set the temperature profile to use when the energy equation
            % is not being solved. The profile must be entered as a pair of arrays
            % specifying positions and temperatures.
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param zFixed:
            %     Array containing positions where the profile is specified
            % :param tFixed:
            %     Array containing temperatures

            ctFunc('flow_setFixedTempProfile', d.domainID, ...
                   zFixed, tFixed);

        end

    end

end

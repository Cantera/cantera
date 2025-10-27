classdef (Abstract) Flow1D < Domain1D
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

    properties (SetAccess = public)
        P  % Flow Pressure. Units: Pa.

        % Boolean flag indicating whether the energy equation is enabled.
        energyEnabled logical

        % Boolean flag indicating whether the diffusive mass fluxes due to the Soret
        % effect is enabled.
        soretEnabled logical

        % ID of the solution object used for calculating transport properties.
        transportModel

        grid  % Grid points from a domain.
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

        %% Domain Properties

        function pressure = get.P(d)
            pressure = ctFunc('flow_pressure', d.domainID);
        end

        function set.P(d, p)
            ctFunc('flow_setPressure', d.domainID, p);
        end

        function flag = get.energyEnabled(d)
            flag = ctFunc('flow_allOfEnergyEnabled', d.domainID);
        end

        function set.energyEnabled(d, flag)
            ctFunc('flow_setEnergyEnabled', d.domainID, flag);
        end

        function set.soretEnabled(d, flag)
            d.soretEnabled = flag;
            ctFunc('flow_enableSoret', d.domainID, flag);
        end

        function model = get.transportModel(d)
            model = ctString('flow_transportModel', d.domainID);
        end

        function set.transportModel(d, model)
            ctFunc('domain_setTransportModel', d.domainID, model);
        end

        function zz = get.grid(d)
            np = d.nPoints;
            zz = ctArray('domain_grid', np, d.domainID);
        end

        function set.grid(d, grid)
            ctFunc('domain_setupGrid', d.domainID, grid);
        end

        %% Flow1D Class Methods

        function setupUniformGrid(d, points, length, start)
            % Set up the solution grid. ::
            %
            %     >> d.setupUniformGrid(points, length, start)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param points:
            %     Number of grid points
            % :param length:
            %     Length of domain
            % :param start:
            %     Start position of domain
            ctFunc('domain_setupUniformGrid', d.domainID, points, length, start);
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

            v = ctArray('domain_values', d.nPoints, d.domainID, component);

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

        function setRefineCriteria(d, ratio, slope, curve, prune)
            % Set the criteria used to refine the grid. ::
            %
            %     >> d.setRefineCriteria(ratio, slope, curve, prune)
            %
            % :param d:
            %    Instance of class :mat:class:`Domain1D`.
            % :param ratio:
            %    Maximum size ratio between adjacent cells.
            % :param slope:
            %    Maximum relative difference in value between adjacent points.
            % :param curve:
            %    Maximum relative difference in slope between adjacent cells.
            % :param prune:
            %    Minimum value for slope or curve for which points will be
            %    retained or curve value is below prune for all components,
            %    it will be deleted, unless either neighboring point is
            %    already marked for deletion.

            if nargin < 3
                ratio = 10.0;
            end

            if nargin < 4
                slope = 0.8;
            end

            if nargin < 5
                curve = 0.8;
            end

            if nargin < 6
                prune = -0.1;
            end

            ctFunc('domain_setRefineCriteria', d.domainID, ratio, slope, curve, prune);
        end

    end

end

classdef (Abstract) Flow < ct.oneD.Domain
    % Flow Class.
    %
    % Base class for objects representing flow domains. The constructor is called
    % by derived classes and cannot be used directly.
    %
    % :param type:
    %    String type of domain. Possible values are:
    %      - `axisymmetric-flow`
    %      - `free-flow`
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String, ID of the flow.

    properties (SetAccess = public)
        P  % Flow Pressure [Pa].

        % Boolean flag indicating whether the energy equation is enabled.
        energyEnabled logical

        % Boolean flag indicating whether the diffusive mass fluxes due to the Soret
        % effect is enabled.
        soretEnabled logical

        % Transport model used for calculating transport properties.
        transportModel string

        grid  % Grid points from a domain.
    end

    methods

        %% Flow Class Constructor

        function obj = Flow(type, phase, name)
            arguments
                type (1,1) string
                phase (1,1) ct.Solution
                name (1,1) string
            end
            id = ct.impl.call('mDomain_newFlow1D', type, phase.solnID, name);
            obj@ct.oneD.Domain(id);
            obj.energyEnabled = false;
            obj.soretEnabled = false;

        end

        %% Domain Properties

        function pressure = get.P(obj)
            pressure = ct.impl.call('mFlow_pressure', obj.domainID);
        end

        function set.P(obj, p)
            ct.impl.call('mFlow_setPressure', obj.domainID, p);
        end

        function flag = get.energyEnabled(obj)
            flag = ct.impl.call('mFlow_allOfEnergyEnabled', obj.domainID);
        end

        function set.energyEnabled(obj, flag)
            ct.impl.call('mFlow_setEnergyEnabled', obj.domainID, flag);
        end

        function set.soretEnabled(obj, flag)
            obj.soretEnabled = flag;
            ct.impl.call('mFlow_enableSoret', obj.domainID, flag);
        end

        function model = get.transportModel(obj)
            model = ct.impl.getString('mFlow_transportModel', obj.domainID);
        end

        function set.transportModel(obj, model)
            ct.impl.call('mDomain_setTransportModel', obj.domainID, model);
        end

        function zz = get.grid(obj)
            np = obj.nPoints;
            zz = ct.impl.getArray('mDomain_grid', np, obj.domainID);
        end

        function set.grid(obj, grid)
            ct.impl.call('mDomain_setupGrid', obj.domainID, grid);
        end

        %% Flow Class Methods

        function setupUniformGrid(obj, points, length, start)
            % Set up the solution grid. ::
            %
            %     >> d.setupUniformGrid(points, length, start)
            %
            % :param points:
            %     Number of grid points
            % :param length:
            %     Length of domain
            % :param start:
            %     Start position of domain
            ct.impl.call('mDomain_setupUniformGrid', obj.domainID, points, length, start);
        end

        function v = values(obj, component)
            % Get the value of a component in a domain. ::
            %
            %     >> d.values(component)
            %
            % :param component:
            %    String component for which the solution is desired.
            % :return:
            %    Value of the component in domain d.

            v = ct.impl.getArray('mDomain_values', obj.nPoints, obj.domainID, component);

        end

        function setProfile(obj, component, zProfile, vProfile)
            % Specify a profile for one component. ::
            %
            %     >> d.setProfile(component, zProfile, vProfile)
            %
            % The solution vector values for this component will be linearly
            % interpolated from the discrete function defined by two input arrays.
            % This method can be called at any time, but is usually used to set the
            % initial guess for the solution.
            %
            % Example (assuming 'd' is an instance of class :mat:class:`ct.oneD.Domain`):
            %    >> zr = [0.0, 0.1, 0.2, 0.4, 0.8, 1.0];
            %
            %    >> v = [500, 650, 700, 730, 800, 900];
            %
            %    >> d.setProfile('T', zr, v);
            %
            % :param component:
            %    Component name.
            % :param zProfile:
            %    Relative (normalized) positions. "zProfile(1) = 0.0" corresponds to the
            %    leftmost grid point in the specified domain, and "zProfile(n) = 1.0"
            %    corresponds to the rightmost grid point
            % :param vProfile:
            %    Array containing component values at positions.

            ct.impl.call('mDomain_setProfile', obj.domainID, ...
                   component, zProfile, vProfile);

        end

        function setFlatProfile(obj, component, v)
            % Set a component to a value across the entire domain. ::
            %
            %     >> d.setFlatProfile(component, v)
            %
            % :param component:
            %    Component to be set.
            % :param v:
            %    Double value to be set.

            ct.impl.call('mDomain_setFlatProfile', obj.domainID, component, v);
        end

        function setFixedTempProfile(obj, zFixed, tFixed)
            % Set a fixed temperature profile. ::
            %
            %     >> d.setFixedTempProfile(zFixed, tFixed)
            %
            % Set the temperature profile to use when the energy equation
            % is not being solved. The profile must be entered as a pair of arrays
            % specifying positions and temperatures.
            %
            % :param zFixed:
            %     Array containing positions where the profile is specified
            % :param tFixed:
            %     Array containing temperatures

            ct.impl.call('mFlow_setFixedTempProfile', obj.domainID, ...
                   zFixed, tFixed);

        end

        function setRefineCriteria(obj, ratio, slope, curve, prune)
            % Set the criteria used to refine the grid. ::
            %
            %     >> d.setRefineCriteria(ratio, slope, curve, prune)
            %
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
            arguments
                obj
                ratio (1,1) double {mustBePositive} = 10.0
                slope (1,1) double {mustBePositive} = 0.8
                curve (1,1) double {mustBePositive} = 0.8
                prune (1,1) double = -0.1
            end

            ct.impl.call('mDomain_setRefineCriteria', obj.domainID, ratio, slope, curve, prune);
        end

    end

end

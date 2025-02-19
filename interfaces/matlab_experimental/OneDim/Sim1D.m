classdef Sim1D < handle
    % Sim1D Class. ::
    %
    %     >> s = Sim1D(domains)
    %
    % A Sim1D object is a container for one-dimensional domains,
    % which are instances of :mat:class:`Domain1D`. The domains are of two
    % types - extended domains, and connector domains.
    %
    % See also: :mat:class:`Domain1D`
    %
    % :param domains:
    %     Cell array of instances of :mat:class:`Domain1D` and its subclasses.
    % :return:
    %     Instance of class :mat:class:`Sim1D`.

    properties (SetAccess = immutable)

        stID % ID of the Sim1D object.

        domains % Domain instances contained within the :mat:class:`Sim1D` object.

    end

    methods

        %% Sim1D Class Constructor

        function s = Sim1D(domains)
            % Create a :mat:class:`Sim1D` object.

            ctIsLoaded;

            s.stID = -1;
            s.domains = domains;

            nd = length(domains);
            ids = zeros(1, nd);

            for n = 1:nd
                ids(n) = domains{n}.domainID;
            end

            s.stID = ctFunc('sim1D_new', nd, ids);

        end

        %% Sim1D Class Destructor

        function delete(s)
            % Delete the :mat:class:`Sim1D` object.

            if isempty(s.stID)
                return
            end
            ctFunc('sim1D_del', s.stID);
        end

        %% Sim1D Utility Methods

        function display(s, fname)
            % Show all domains.

            if nargin == 1
                fname = '-';
            end

            ctFunc('sim1D_show', s.stID, fname);
        end

        function restore(s, fname, id)
            % Restore a previously-saved solution. ::
            %
            %     >> s.restore(fname, id)
            %
            % This method can be used to provide an initial guess for the solution.
            %
            % See also: :mat:class:`save`
            %
            % :param s:
            %     Instance of class :mat:class:`Sim1D`.
            % :param fname:
            %     File name of an YAML or HDF file containing solution information.
            % :param id:
            %     ID of the element that should be restored.

            ctFunc('sim1D_restore', s.stID, fname, id)
        end

        function save(s, fname, id, desc)
            % Save a solution to a file. ::
            %
            %     >> s.save(fname, id, desc)
            %
            % The output file is in a format that can be used by :mat:class:`restore`
            %
            % :param s:
            %     Instance of class :mat:class:`Sim1D`.
            % :param fname:
            %     File name where YAML or HDL file should be written.
            % :param id:
            %     ID to be assigned to the file element when it is written.
            % :param desc:
            %     Description to be written to the output file.

            if nargin < 3
                error('Not enough input arguments');
            elseif nargin == 3
                desc = '';
            end

            ctFunc('sim1D_save', s.stID, fname, id, desc);
        end

        function x = getSolution(s, domain, component)
            % Get a solution component in one domain. ::
            %
            %     >> s.getSolution(domain, component)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param domain:
            %    String name of the domain from which the solution is desired.
            % :param component:
            %    String component for which the solution is desired. If omitted,
            %    solution for all of the components will be returned in
            %    an :math:`nPoints \times nComponents` array.
            % :return:
            %    Either an :math:`nPoints \times 1` vector, or
            %    :math:`nPoints \times nComponents` array.

            idom = s.stackIndex(domain);
            d = s.domains{idom};
            np = d.nPoints;

            if nargin == 3
                icomp = d.componentIndex(component);
                x = zeros(1, np);

                for n = 1:np
                    x(n) = ctFunc('sim1D_value', s.stID, idom - 1, icomp - 1, n - 1);
                end

            else
                nc = d.nComponents;
                x = zeros(nc, np);

                for m = 1:nc

                    for n = 1:np
                        x(m, n) = ctFunc('sim1D_value', s.stID, idom - 1, m - 1, n - 1);
                    end

                end

            end

        end

        function solve(s, loglevel, refineGrid)
            % Solve the problem. ::
            %
            %     >> s.solve(loglevel, refineGrid)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param loglevel:
            %    Integer flag controlling the amount of diagnostic output.
            %    Zero suppresses all output, and 5 produces very verbose output.
            % :param refineGrid:
            %    Integer, 1 to allow grid refinement, 0 to disallow.

            ctFunc('sim1D_solve', s.stID, loglevel, refineGrid);
        end

        function writeStats(s)
            % Print statistics for the current solution. ::
            %
            %     >> s.writeStats
            %
            % Prints a summary of the number of function and Jacobian evaluations
            % for each grid, and the CPU time spent on each one.
            %
            % :param s:
            %     Instance of class class :mat:class:`Sim1D`

            ctFunc('sim1D_writeStats', s.stID, 1);
        end

        %% Sim1D Get Methods

        function getInitialSoln(s)
            % Get the initial solution.

            ctFunc('sim1D_getInitialSoln', s.stID);
        end

        function n = stackIndex(s, name)
            % The index of a domain in a Sim1D given its name. ::
            %
            %     >> n = s.stackIndex(name)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param name:
            %    If double, the value is :returned. Otherwise, the name is
            %    looked up and its index is :returned.
            % :return:
            %    Index of domain.

            if isa(name, 'double')
                n = name;
            else
                n = ctFunc('sim1D_domainIndex', s.stID, name);

                if n >= 0
                    n = n + 1;
                else
                    error('Domain not found');
                end

            end

        end

        function z = grid(s, name)
            % Get the grid in one domain. ::
            %
            %     >> z = s.grid(name)
            %
            % :param s:
            %     Instance of class :mat:class:`Sim1D`.
            % :param name:
            %     Name of the domain for which the grid should be retrieved.
            % :return:
            %     The grid in domain name.

            n = s.stackIndex(name);
            d = s.domains{n};
            z = d.gridPoints;
        end

        function r = residual(s, domain, rdt, count)
            % Evaluate the multi-domain residual function. ::
            %
            %     >> r = s.residual(domain, rdt, count)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param domain:
            %    Name of the domain.
            % :param rdt:
            %    Reciprocal of the time step. If omitted, the default value is used.
            % :param count:
            %    Set to zero to omit this call from the statistics.
            % :return:
            %    The multi-domain residual function.

            if nargin == 2
                rdt = 0.0;
                count = 0;
            end

            idom = s.stackIndex(domain);
            d = s.domains{idom};

            nc = d.nComponents;
            np = d.nPoints;

            r = zeros(nc, np);
            ctFunc('sim1D_eval', s.stID, rdt, count);

            for m = 1:nc

                for n = 1:np
                    r(m, n) = ctFunc('sim1D_workValue', s.stID, idom - 1, m - 1, n - 1);
                end

            end

        end

        %% Sim1D Set Methods

        function setFixedTemperature(s, T)
            % Set the temperature used to fix the spatial location of a
            % freely propagating flame. ::
            %
            %     >> s.setFixedTemperature(T)
            %
            % :param T:
            %    Double Temperature to be set. Unit: K.

            if T <= 0
                error('temperature must be positive');
            end

            ctFunc('sim1D_setFixedTemperature', s.stID, T);
        end

        function setFlatProfile(s, domain, comp, v)
            % Set a component to a value across the entire domain. ::
            %
            %     >> s.setFlatProfile(domain, comp, v)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param domain:
            %    Integer ID of the domain.
            % :param comp:
            %    Component to be set.
            % :param v:
            %    Double value to be set.

            ctFunc('sim1D_setFlatProfile', s.stID, domain - 1, comp - 1, v);
        end

        function setGridMin(s, domain, gridmin)
            % Set the minimum grid spacing on domain. ::
            %
            %     >> s.setGridMin(domain, gridmin)
            %
            % :param domain:
            %    Integer ID of the domain.
            % :param gridmin:
            %    Double minimum grid spacing.

            ctFunc('sim1D_setGridMin', s.stID, domain - 1, gridmin);
        end

        function setMaxJacAge(s, ss_age, ts_age)
            % Set the number of times the Jacobian will be used before it
            % is recomputed. ::
            %
            %     >> s.setMaxJacAge(ss_age, ts_age)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param ss_age:
            %    Maximum age of the Jacobian for steady state analysis.
            % :param ts_age:
            %    Maximum age of the Jacobian for transient analysis. If not
            %    specified, defaults to 'ss_age'.

            if nargin == 2
                ts_age = ss_age;
            end

            ctFunc('sim1D_setMaxJacAge', s.stID, ss_age, ts_age);
        end

        function setProfile(s, name, comp, p)
            % Specify a profile for one component. ::
            %
            %     >> s.setProfile(name, comp, p)
            %
            % The solution vector values for this component will be
            % linearly interpolated from the discrete function defined by
            % p(:, 1) vs p(:, 2).
            % Note that "p(1, 1) = 0.0" corresponds to the leftmost grid
            % point in the specified domain, and "p(1, n) = 1.0"
            % corresponds to the rightmost grid point. This method can be
            % called at any time, but is usually used to set the initial
            % guess for the solution.
            %
            % Example (assuming 's' is an instance of class :mat:class:`Sim1D`):
            %    >> zr = [0.0, 0.1, 0.2, 0.4, 0.8, 1.0];
            %
            %    >> v = [500, 650, 700, 730, 800, 900];
            %
            %    >> s.setProfile(1, 2, [zr, v]);
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param name:
            %    Domain name.
            % :param comp:
            %    Component number.
            % :param p:
            %    n x 2 array, whose columns are the relative (normalized)
            %    positions and the component values at those points. The
            %    number of positions 'n' is arbitrary.

            if isa(name, 'double')
                n = name;
            else
                n = s.domainIndex(name);
            end

            d = s.domains{n};

            if isa(comp, 'double') || isa(comp, 'cell')
                c = comp;
            elseif isa(comp, 'char')
                c = {comp};
            else
                error('Wrong type.');
            end

            np = length(c);
            sz = size(p);

            if sz(1) == np + 1;

                for j = 1:np
                    ic = d.componentIndex(c{j});
                    ctFunc('sim1D_setProfile', s.stID, ...
                            n - 1, ic - 1, sz(2), p(1, :), sz(2), p(j + 1, :));
                end

            elseif sz(2) == np + 1;
                ic = d.componentIndex(c{j});
                ctFunc('sim1D_setProfile', s.stID, ...
                        n - 1, ic - 1, sz(1), p(:, 1), sz(1), p(:, j + 1));
            else
                error('Wrong profile shape.');
            end

        end

        function setRefineCriteria(s, n, ratio, slope, curve, prune)
            % Set the criteria used to refine the grid. ::
            %
            %     >> s.setRefineCriteria(n, ratio, slope, curve, prune)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
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

            ctFunc('sim1D_setRefineCriteria', s.stID, ...
                    n - 1, ratio, slope, curve, prune);
        end

        function setTimeStep(s, stepsize, steps)
            % Specify a sequence of time steps. ::
            %
            %     >> s.setTimeStep(stepsize, steps)
            %
            % :param stepsize:
            %    Initial step size.
            % :param steps:
            %    Vector of number of steps to take before re-attempting solution
            %    of steady-state problem.
            %    For example, steps = [1, 2, 5, 10] would cause one timestep to be
            %    taken first time the steady-state solution attempted.
            %    If this failed, two time steps would be taken.

            ctFunc('sim1D_setTimeStep', s.stID, stepsize, length(steps), steps);
        end

        function setValue(s, n, comp, localPoints, v)
            % Set the value of a single entry in the solution vector. ::
            %
            %     >> s.setValue(n, comp, localPoints, v)
            %
            % Example (assuming 's' is an instance of class :mat:class:`Sim1D`) ::
            %
            %    >> s.setValue(3, 5, 1, 5.6);
            %
            % This sets component 5 at the leftmost point (local point 1)
            % in domain 3 to the value 5.6. Note that the local index
            % always begins at 1 at the left of each domain, independent of
            % the global index of the point, which depends on the location
            % of this domain in the class :mat:class:`Sim1D` object.
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param n:
            %    Domain number.
            % :param comp:
            %    Component number.
            % :param localPoints:
            %    Local index of the grid point in the domain.
            % :param v:
            %    Value to be set.

            ctFunc('sim1D_setValue', s.stID, n - 1, comp - 1, localPoints - 1, v);
        end

    end

end

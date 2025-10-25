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

            s.stID = ctFunc('sim1D_newSim1D', ids);

        end

        %% Sim1D Class Destructor

        function delete(s)
            % Delete the :mat:class:`Sim1D` object.
            ctFunc('sim1D_del', s.stID);
        end

        %% Sim1D Utility Methods

        function show(s)
            % Show all domains.
            ctFunc('sim1D_show', s.stID);
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

        function n = domainIndex(s, name)
            % The index of a domain in a Sim1D given its name. ::
            %
            %     >> n = s.domainIndex(name)
            %
            % :param s:
            %    Instance of class :mat:class:`Sim1D`.
            % :param name:
            %    If double, the value is :returned. Otherwise, the name is
            %    looked up and its index is :returned.
            % :return:
            %    Index of domain.

            n = ctFunc('sim1D_domainIndex', s.stID, name);

            if n >= 0
                n = n + 1;
            else
                error('Domain not found');
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

            ctFunc('sim1D_setTimeStep', s.stID, stepsize, steps);
        end

    end

end

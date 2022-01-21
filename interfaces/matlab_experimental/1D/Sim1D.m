classdef Stack < handle

    properties
        st_id
        domains
    end

    methods
        %% Stack class constructor

        function s = Stack(domains)
            % A stack object is a container for one-dimensional domains,
            % which are instances of class Domain1D. The domains are of two
            % types - extended domains, and connector domains.
            %
            % parameter domains:
            %    Vector of domain instances.
            % return:
            %    Instance of class 'Stack'.

            checklib;

            s.st_id = -1;
            s.domains = domains;
            if nargin == 1
                nd = length(domains);
                ids = zeros(1, nd);
                for n=1:nd
                    ids(n) = domains(n).dom_id;
                end
                s.st_id = calllib(ct, 'sim1D_new', nd, ids);
            else
                help(Stack);
                error('Wrong number of parameters.');
            end
%             if s.st_id < 0
%                 error(geterr);
%             end
        end

        %% Utility Methods

        function st_clear(s)
            % Delete the Sim1D object
            checklib;
            calllib(ct, 'sim1D_del', s.st_id);
        end

        function display(s, fname)
            % Show all domains.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter fname:
            %    File to write summary to. If omitted, output is to the
            %    command window.

            checklib;
            if nargin == 1
                fname = '-';
            end
            calllib(ct, 'sim1D_showSolution', s.st_id, fname);
        end

        %% Stack Methods

        function n = stackIndex(s, name)
            % Get the index of a domain in a stack given its name.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter name:
            %    If double, the value is returned. Otherwise, the name is
            %    looked up and its index is returned.
            % return:
            %    Index of domain.

            checklib;
            if isa(name, 'double')
                n = name;
            else
                n = calllib(ct, 'sim1D_domainIndex', s.st_id, name);
                if n >= 0
                    n = n+1;
                else
                    error('Domain not found');
                end
            end
        end

        function z = grid(s, name)
            % Get the grid in one domain.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter name:
            %    Name of the domain for which the grid should be retrieved.
            % return:
            %    The grid in domain name.

            n = s.stackIndex(name);
            d = s.domains(n);
            z = d.gridPoints;
        end

        function plotSolution(s, domain, component)
            % Plot a specified solution component.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter domain:
            %    Name of domain from which the component should be
            %    retrieved.
            % parameter component:
            %    Name of the component to be plotted.

            n = s.stackIndex(domain);
            d = s.domains(n);
            z = d.gridPoints;
            x = s.solution(domain, component);
            plot(z, x);
            xlabel('z (m)');
            ylabel(component);
        end

        function r = resid(s, domain, rdt, count)
            % Get the residuals.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter domain:
            %    Name of the domain.
            % parameter rdt:
            % parameter count:
            % return:

            checklib;

            if nargin == 2
                rdt = 0.0;
                count = 0;
            end

            idom = s.stackIndex(domain);
            d = s.domains(idom);

            nc = d.nComponents;
            np = d.nPoints;

            r = zeros(nc, np);
            calllib(ct, 'sim1D_eval', s.st_id, rdt, count);
            for m = 1:nc
                for n = 1:np
                    r(m, n) = calllib(ct, 'sim1D_workValue', ...
                                      s.st_id, idom - 1, m - 1, n - 1);
                end
            end
        end

        function restore(s, fname, id)
            % Restore a previously-saved solution.
            % This method can be used ot provide an initial guess for the
            % solution.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter fname:
            %    File name of an XML file containing solution info.
            % parameter id:
            %    ID of the element that should be restored.
            checklib;
            calllib(ct, 'sim1D_restore', s.st_id, fname, id)
        end

        function saveSoln(s, fname, id, desc)
            % Save a solution to a file.
            % The output file is in a format that can be used by 'restore'.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter fname:
            %    File name where XML file should be written.
            % parameter id:
            %    ID to be assigned to the XMl element when it is written.
            % parameter desc:
            %    Description to be written to the output file.
            checklib;

            if nargin == 1
                fname = 'soln.xml';
                id = 'solution';
                desc = '--';
            elseif nargin == 2
                id = 'solution';
                desc = '--';
            elseif nargin == 3
                desc = '--';
            end
            calllib(ct, 'sim1D_save', s.st_id, fname, id, desc);
        end

        function setFlatProfile(s, domain, comp, v)
            % Set a component to a value across the entire domain.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter domain:
            %    Integer ID of the domain.
            % parameter comp:
            %    Component to be set.
            % parameter v:
            %    Double value to be set.
            checklib;
            calllib(ct, 'sim1D_setFlatProfile', s.st_id, ...
                    domain - 1, comp - 1, v);
        end

        function setMaxJacAge(s, ss_age, ts_age)
            % Set the number of times the Jacobian will be used before it
            % is recomputed.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter ss_age:
            %    Maximum age of the Jacobian for steady state analysis.
            % parameter ts_age:
            %    Maximum age of the Jacobian for transient analysis. If not
            %    specified, deftauls to 'ss_age'.
            checklib;

            if nargin == 2
                ts_age = ss_age;
            end
            calllib(ct, 'sim1D_setMaxJacAge', s.st_id, ss_age, ts_age);
        end

        function setProfile(s, name, comp, p)
            % Specify a profile for one component,
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
            % Example (assuming 's' is an instance of class 'Stack'):
            %    >> zr = [0.0, 0.1, 0.2, 0.4, 0.8, 1.0];
            %    >> v = [500, 650, 700, 730, 800, 900];
            %    >> s.setProfile(1, 2, [zr, v]);
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter name:
            %    Domain name.
            % parameter comp:
            %    Component number.
            % parameter p:
            %    n x 2 array, whose columns are the relative (normalized)
            %    positions and the component values at those points. The
            %    number of positions 'n' is arbitrary.

            checklib;

            if isa(name, 'double')
                n = name;
            else
                n = s.domainIndex(name);
            end

            d = s.domains(n);

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
                    calllib(ct, 'sim1D_setProfile', s.st_id, ...
                            n - 1, ic - 1, sz(1), p(1, :), sz(1), p(j+1, :));
                end
            elseif sz(2) == np + 1;
                ic = d.componentIndex(c{j});
                calllib(ct, 'sim1D_setProfile', s.st_id, ...
                        n - 1, ic - 1, sz(2), p(:, 1), sz(2), p(:, j+1));
            else
                error('Wrong profile shape.');
            end
        end

        function setRefineCriteria(s, n, ratio, slope, curve, prune)
            % Set the criteria used to refine the grid.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter ratio:
            %    Maximum size ratio between adjacent cells.
            % parameter slope:
            %    Maximum relative difference in value between adjacent
            %    points.
            % parameter curve:
            %    Maximum relative difference in slope between adjacent
            %    cells.
            % parameter prune:
            %    Minimum value for slope or curve for which points will be
            %    retained or curve value is below prune for all components,
            %    it will be deleted, unless either neighboring point is
            %    already marked for deletion.

            checklib;

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
            calllib(ct, 'sim1D_setRefineCriteria', s.st_id, ...
                    n - 1, ratio, slope, curve, prune);
        end

        function setTimeStep(s, stepsize, steps)
            % Specify a sequence of time steps.
            %
            % parameter stepsize:
            %    Initial step size.
            % parameter steps:
            %    Vector of number of steps to take before re-attempting
            %    solution of steady-state problem.
            %    For example, steps = [1, 2, 5, 10] would cause one time
            %    step to be taken first time the steady-state solution
            %    attempted. If this failed, two time steps would be taken.
            checklib;
            calllib(ct, 'sim1D_', s.st_id, ...
                    stepsize, length(steps), steps);
        end

        function setValue(s, n, comp, localPoints, v)
            % Set the value of a single entry in the solution vector.
            %
            % Example (assuming 's' is an instance of class 'Stack'):
            %
            %    setValue(s, 3, 5, 1, 5.6);
            %
            % This sets component 5 at the leftmost point (local point 1)
            % in domain 3 to the value 5.6. Note that the local index
            % always begins at 1 at the left of each domain, independent of
            % the global index of the point, wchih depends on the location
            % of this domain in the stack.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter n:
            %    Domain number.
            % parameter comp:
            %    Component number.
            % parameter localPoints:
            %    Local index of the grid point in the domain.
            % parameter v:
            %    Value to be set.
            checklib;
            calllib(ct, 'sim1D_setValue', s.st_id, ...
                    n - 1, comp -  1, localPoints - 1, v);
        end

        function x = solution(s, domain, component)
            % Get a solution component in one domain.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter domain:
            %    String name of the domain from which the solution is
            %    desired.
            % parameter component:
            %    String component for which the solution is desired. If
            %    omitted, solution for all of the components will be
            %    returned in an 'nPoints' x 'nComponnts' array.
            % return:
            %    Either an 'nPoints' x 1 vector, or 'nPoints' x
            %    'nCOmponents' array.
            checklib;

            idom = s.stackIndex(domain);
            d = s.domains(idom);
            np = d.nPoints;
            if nargin == 3
                icomp = d.componentIndex(component);
                x = zeros(1, np);
                for n = 1:np
                    x(n) = calllib(ct, 'sim1D_value', s.st_id, ...
                                   idom - 1, icomp - 1, n - 1);
                end
            else
                nc = d.nComponents;
                x = zeros(nc, np);
                for m = 1:nc
                    for n = 1:np
                        x(m, n) = calllib(ct, 'sim1D_value', s.st_id, ...
                                          idom - 1, m - 1, n - 1);
                    end
                end
            end
        end

        function solve(s, loglevel, refine_grid)
            % Solve the problem.
            %
            % parameter s:
            %    Instance of class 'Stack'.
            % parameter loglevel:
            %    Integer flag controlling the amount of diagnostic output.
            %    Zero supresses all output, and 5 produces very verbose
            %    output.
            % parameter refine_grid:
            %    Integer, 1 to allow grid refinement, 0 to disallow.
            checklib;
            calllib(ct, 'sim1D_solve', s.st_id, loglevel, refine_grid);
        end

%         function b = subsref(s, index)
%             % Redefine subscripted references.
%             switch index.type
%                 case '()'
%                     b = s.domains(index.subs{:});
%                 case '.'
%                     n = s.domainIndex(index.subs);
%                     b = s.domains(n);
%                 otherwise
%                     error('syntax error');
%             end
%         end

        function writeStats(s)
            % Print statistics for the current solution.
            % Prints a summary of the number of function and Jacobian
            % evaluations for each grid, and the CPU time spent on each
            % one.
            checklib;
            calllib(ct, 'sim1D_writeStats', s.st_id, 1);
        end

    end
end

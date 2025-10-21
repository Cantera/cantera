classdef Interface < Solution
    % Interface Class ::
    %
    %     >> s = Interface(src, name, p1, p2)
    %
    % See :ref:`sec-yaml-ideal-surface` and :ref:`sec-yaml-guide-adjacent`.
    %
    % :param src: YAML file containing the interface or edge phase.
    % :param name: Name of the interface or edge phase in the YAML file.
    % :param varargin:
    %     Optional list of phases pi adjacent to the interface; if omitted, adjacent
    %     phases are added automatically.
    % :return:
    %     Instance of class :mat:class:`Interface`.

    properties (SetAccess = public)

        % Surface coverages of the species on an interface.
        % Unit: kmol/m^2 for surface phases, kmol/m for edge phases.
        siteDensity

        coverages % Surface coverages of the species on an interface.

    end

    properties (SetAccess = protected)
        concentrations % Concentrations of the species on an interface.
        nAdjacent % Number of adjacent phases.
        adjacentNames % Names of adjacent phases.
    end

    methods
        %% Interface Class Constructor

        function s = Interface(src, name, varargin)
            % Create an :mat:class:`Interface` object.

            ctIsLoaded;

            na = nargin - 2;

            % Get ID of adjacent phases
            adj = [];
            for i = 1:na
                adj(i) = varargin{i}.solnID;
            end

            ID = ctFunc('sol_newInterface', src, name, adj);

            % Inherit methods and properties from Solution
            s@Solution(ID);
            s.nAdjacent = ctFunc('sol_nAdjacent', ID);
            s.adjacentNames = {};
            for i = 1:s.nAdjacent
                s.adjacentNames{i} = ctString('sol_adjacentName', ID, i-1);
            end
        end

        %% Interface Class Destructor

        function delete(s)
            % Delete :mat:class:`Interface` object.

            delete@Solution(s);
        end

        %% Interface Get Methods

        function adj = adjacent(s, name)
            % Get adjacent phase of an interface by name.
            exact_match = strcmp(s.adjacentNames, name);
            if sum(exact_match) ~= 1
                error(['No adjacent phase with name ''' name ''' found.'])
            end
            location = find(exact_match);
            id = ctFunc('sol_adjacent', s.solnID, location-1);
            adj = Solution(id);
        end

        function c = get.coverages(s)
            nsp = s.nSpecies;
            c = ctArray('surf_getCoverages', nsp, s.tpID);
        end

        function d = get.siteDensity(s)
            d = ctFunc('surf_siteDensity', s.tpID);
        end

        function c = get.concentrations(s)
            nsp = s.nSpecies;
            c = ctArray('thermo_getConcentrations', nsp, s.tpID);
        end

        function set.coverages(s, cov)
            % Set the fraction of sites covered by each species.
            %
            % s.coverages = cov
            %
            % :param s:
            %      Instance of class :mat:class:`Interface`
            % :param cov:
            %      Vector or string coverage of the species.

            nsp = s.nSpecies;
            if isa(cov, 'double')
                if length(cov) ~= nsp
                    error('wrong size for coverage array');
                end

                ctFunc('surf_setCoverages', s.tpID, cov);
            elseif isa(cov, 'char')
                ctFunc('surf_setCoveragesByName', s.tpID, cov);
            end

        end

        function set.siteDensity(s, d)
            % Set the site density of a phase on an interface.
            %
            % s.siteDensity = d
            %
            % :param s:
            %      Instance of class :mat:class:`Interface`
            % :param d
            %    Double site density. Unit: kmol/m^2 for surface phases,
            %    kmol/m for edge phases.
            ctFunc('surf_setSiteDensity', s.tpID, d);
        end

        function setUnnormalizedCoverages(s, cov)
            % Set surface coverages without normalizing to force sumo(cov) == 1.0.
            % This should be used only when calculating partial derivatives
            % with respect to cov[k] by finite difference.
            %
            % s.setUnnormalizedCoverages(cov)
            %
            % :param s:
            %      Instance of class :mat:class:`Interface`
            % :param cov:
            %      Vector coverage of the species.

            nsp = s.nSpecies;
            if isa(cov, 'double')
                if length(cov) ~= nsp
                    error('wrong size for coverage array');
                end

                ctFunc('surf_setCoverages', s.tpID, cov, 0);
            else
                error('Coverage must be a numeric array');
            end

        end

    end

end

classdef Interface < ct.Solution
    % Interface Class ::
    %
    %     >> s = ct.Interface(src, name, p1, p2)
    %
    % See :ref:`sec-yaml-ideal-surface` and :ref:`sec-yaml-guide-adjacent`.
    %
    % :param src: YAML file containing the interface or edge phase.
    % :param name: Name of the interface or edge phase in the YAML file.
    % :param varargin:
    %     Optional list of phases pi adjacent to the interface; if omitted, adjacent
    %     phases are added automatically.

    properties (SetAccess = public)

        % Density of sites for interface species [kmol/m² for surface phases or
        % kmol/m for edge phases].
        siteDensity

        coverages % Surface coverages of the species on an interface.

    end

    properties (SetAccess = protected)
        nAdjacent % Number of adjacent phases.
        adjacentNames % Names of adjacent phases.
    end

    methods
        %% Interface Class Constructor

        function obj = Interface(src, name, varargin)
            ct.isLoaded(true);

            na = nargin - 2;

            % Get ID of adjacent phases
            adj = [];
            for i = 1:na
                adj(i) = varargin{i}.solnID;
            end

            ID = ct.impl.call('mSol_newInterface', src, name, adj);

            % Inherit methods and properties from Solution
            obj@ct.Solution(ID);
            obj.nAdjacent = ct.impl.call('mSol_nAdjacent', ID);
            obj.adjacentNames = {};
            for i = 1:obj.nAdjacent
                obj.adjacentNames{i} = ct.impl.getString('mSol_adjacentName', ID, i-1);
            end
        end

        %% Interface Get Methods

        function adj = adjacent(obj, name)
            % Get adjacent phase of an interface by name.
            id = ct.impl.call('mSol_adjacentByName', obj.solnID, name);
            adj = ct.Solution(id);
        end

        function c = get.coverages(obj)
            nsp = obj.nSpecies;
            c = ct.impl.getArray('mSurf_getCoverages', nsp, obj.tpID);
        end

        function d = get.siteDensity(obj)
            d = ct.impl.call('mSurf_siteDensity', obj.tpID);
        end

        function set.coverages(obj, cov)
            nsp = obj.nSpecies;
            if isa(cov, 'double')
                if length(cov) ~= nsp
                    error('wrong size for coverage array');
                end

                ct.impl.call('mSurf_setCoverages', obj.tpID, cov);
            elseif isa(cov, 'char')
                ct.impl.call('mSurf_setCoveragesByName', obj.tpID, cov);
            end

        end

        function set.siteDensity(obj, d)
            ct.impl.call('mSurf_setSiteDensity', obj.tpID, d);
        end

        function setUnnormalizedCoverages(obj, cov)
            % Set surface coverages without normalizing to force ``sum(cov) == 1.0``.
            % This should be used only when calculating partial derivatives
            % with respect to cov[k] by finite difference.
            %
            % s.setUnnormalizedCoverages(cov)
            %
            % :param cov:
            %      Vector coverage of the species.

            nsp = obj.nSpecies;
            if isa(cov, 'double')
                if length(cov) ~= nsp
                    error('wrong size for coverage array');
                end

                ct.impl.call('mSurf_setCoverages', obj.tpID, cov, 0);
            else
                error('Coverage must be a numeric array');
            end

        end

    end

end

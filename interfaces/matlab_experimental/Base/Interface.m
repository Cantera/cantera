classdef Interface < handle & ThermoPhase & Kinetics
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

    properties (SetAccess = immutable)
        solnID % ID of the interface.
        interfaceName % Name of the interface.
    end

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

            ID = ctFunc('soln_newInterface', src, name, na, adj);

            % Inherit methods and properties from ThermoPhase and Kinetics
            s@ThermoPhase(ID);
            s@Kinetics(ID);
            s.solnID = ID;
            s.interfaceName = name;
            s.nAdjacent = ctFunc('soln_nAdjacent', ID);
            s.adjacentNames = {};
            for i = 1:s.nAdjacent
                s.adjacentNames{i} = ctString('soln_adjacentName', ID, i-1);
            end
        end

        %% Interface Class Destructor

        function delete(s)
            % Delete :mat:class:`Interface` object.
            ctFunc('soln_del', s.solnID);
        end

        %% Interface Get Methods

        function adj = adjacent(s, name)
            % Get adjacent phase of an interface by name.
            exact_match = strcmp(s.adjacentNames, name);
            if sum(exact_match) ~= 1
                error(['No adjacent phase with name ''' name ''' found.'])
            end
            location = find(exact_match);
            id = ctFunc('soln_adjacent', s.solnID, location-1);
            adj = Solution(id);
        end

        function c = get.coverages(s)
            surfID = s.tpID;
            nsp = s.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('surf_getCoverages', surfID, pt);
            c = pt.Value;
        end

        function d = get.siteDensity(s)
            surfID = s.tpID;
            d = ctFunc('surf_siteDensity', surfID);
        end

        function c = get.concentrations(s)
            surfID = s.tr_id;
            nsp = s.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('surf_getConcentrations', surfID, xx);
            c = pt.Value;
        end

        function set.coverages(s, val)
            % Set surface coverages of the species on an interface.
            %
            % s.coverages = {cov, norm}
            %
            % :param s:
            %      Instance of class :mat:class:`Interface`
            % :param cov:
            %      Coverage of the species. ``cov`` can be either a vector of
            %      length ``n_surf_species``, or a string in the format
            %      ``'Species:Coverage, Species:Coverage'``
            % :param norm:
            %      Optional flag that denotes whether or not to normalize the species
            %      coverages. ``norm`` is either of the two strings ``'nonorm'``` or
            %      ``'norm'``. If left unset, the default is `norm`. This only works if
            %      ``s`` is a vector, not a string. Since unnormalized coverages can lead to
            %      unphysical results, ``'nonorm'`` should be used only in rare cases, such
            %      as computing partial derivatives with respect to a species coverage.

            if iscell(val) && numel(val) >= 1 && numel(val) <= 2
                cov = val{1};
                norm_flag = 1;

                if numel(val) == 2
                    norm = val{2};
                    if strcmp(norm, 'nonorm')
                        norm_flag = 0;
                    end
                end
            else
                error('Input must be a cell array {cov} or {cov, norm}');
            end

            surfID = s.tpID;
            nsp = s.nSpecies;
            [m, n] = size(cov);

            if isa(cov, 'double')
                sz = length(cov);

                if sz ~= nsp
                    error('wrong size for coverage array');
                end

                if ((m == nsp && n == 1) || (m == 1 & n == nsp))
                    ctFunc('surf_setCoverages', surfID, cov, norm_flag);
                else error('wrong size for coverage array');
                end

            elseif isa(cov, 'char')
                ctFunc('surf_setCoveragesByName', surfID, cov);
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

            surfID = s.tpID;
            ctFunc('surf_setSiteDensity', surfID, d);
        end

    end

end

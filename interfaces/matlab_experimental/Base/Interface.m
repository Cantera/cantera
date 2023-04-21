classdef Interface < handle & ThermoPhase & Kinetics
    % Interface Class ::
    %
    %     >> s = Interface(src, id, p1, p2, p3, p4)
    %
    % See `ideal-surface <https://cantera.org/documentation/docs-2.6/sphinx/html/yaml/phases.html#sec-yaml-ideal-surface>`__
    % and `Declaring adjacent phases <https://cantera.org/tutorials/yaml/phases.html#declaring-adjacent-phases>`__.
    %
    % :param varagin:
    %     Variable number of inputs consisting of the following:
    %       - src: YAML file containing the interface or edge phase.
    %       - id: Name of the interface or edge phase in the YAML file.
    %     Optional:
    %       - p1: 1st Adjoining phase to the interface.
    %       - p2: 2nd Adjoining phase to the interface.
    %       - p3: 3rd Adjoining phase to the interface.
    %       - p4: 4th Adjoining phase to the interface.
    % :return:
    %     Instance of class :mat:class:`Interface`.

    properties (SetAccess = public)

        % Surface coverages of the species on an interface.
        % Unit: kmol/m^2 for surface phases, kmol/m for edge phases.
        siteDensity

    end

    properties (SetAccess = protected)
        concentrations % Concentrations of the species on an interface.
    end

    methods
        %% Interface Class Constructor

        function s = Interface(varargin)
            % Create an :mat:class:`Interface` object.

            ctIsLoaded;

            src = varargin{1};
            id = varargin{2};

            t = ThermoPhase(src, id);
            s@ThermoPhase(src, id);

            args = varargin(3:end);

            s@Kinetics(t, src, id, args{:});
            s.tpID = t.tpID;
        end

        %% Interface Class Destructor

        function delete(s)
            % Delete :mat:class:`Interface` object.
            disp('Interface class object has been deleted');
        end

        %% Interface Get Methods

        function c = coverages(s)
            % Surface coverages of the species on an interface.

            surfID = s.tpID;
            nsp = s.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('surf_getCoverages', surfID, pt);
            c = pt.Value;
        end

        function d = get.siteDensity(s)
            surfID = s.tpID;
            d = calllibt(ct, 'surf_siteDensity', surfID);
        end

        function c = get.concentrations(s)
            surfID = s.tr_id;
            nsp = s.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('surf_getConcentrations', surfID, xx);
            c = pt.Value;
        end

        function setCoverages(s, cov, norm)
            % Set surface coverages of the species on an interface.
            %
            % s.setCoverages(cov, norm)
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

            if nargin == 3 && strcmp(norm, 'nonorm')
                norm_flag = 0;
            else
                norm_flag = 1;
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

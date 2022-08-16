classdef Interface < handle & ThermoPhase & Kinetics

    properties
        coverages
        siteDensity
    end

    methods
        %% Interface class constructor

        function s = Interface(src, id, p1, p2, p3, p4)
            % INTERFACE  Interface class constructor.
            % s = Interface(src, id, p1, p2, p3, p4)
            %
            % See also: :mat:func:`importEdge`, :mat:func:`importInterface`
            %
            % :param src:
            %     YAML file containing the interface or edge phase.
            % :param id:
            %     Name of the interface or edge phase in the YAML file.
            % :param p1:
            %     Adjoining phase to the interface.
            % :param p2:
            %     Adjoining phase to the interface.
            % :param p3:
            %     Adjoining phase to the interface.
            % :param p4:
            %     Adjoining phase to the interface.
            % :return:
            %     Instance of class :mat:func:`Interface`
            %
            checklib;
            t = ThermoPhase(src, id);
            s@ThermoPhase(src, id);
            if nargin == 2
                args = {};
            elseif nargin == 3
                args = {p1};
            elseif nargin == 4
                args = {p1, p2};
            elseif nargin == 5
                args = {p1, p2, p3};
            elseif nargin == 6
                args = {p1, p2, p3, p4};
            end
            s@Kinetics(t, src, id, args{:});
            s.tpID = t.tpID;
        end

        %% Interface methods

        function c = get.coverages(s)
            % GET.COVERAGES  Get the surface coverages of the species on an interface.
            % c = s.coverages
            % :param s:
            %     Instance of class :mat:func:`Interface` with surface species
            % :return:
            %     If no output value is assigned, a bar graph will be plotted.
            %     Otherwise, a vector of length ``n_surf_species`` will be
            %     returned.
            %
            surfID = s.tpID;
            nsp = s.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('surf_getCoverages', surfID, xx);
            c = pt.Value;
        end

        function d = get.siteDensity(s)
            % GET.SITEDENSITY  Get the surface coverages of the species on an interface.
            % c = s.siteDensity
            % :param s:
            %     Instance of class :mat:func:`Interface` with surface species
            % :return:
            %    Double site density. Unit: kmol/m^2 for surface phases,
            %    kmol/m for edge phases.
            %
            surfID = s.tpID;
            d = calllibt(ct, 'surf_siteDensity', surfID);
        end

        function c = concentrations(s)
            % Get the concentrations of the species on an interface.
            %
            % :return:
            %    Vector of length "n_surf_species" for concentration.

            surfID = s.tr_id;
            nsp = s.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('surf_getConcentrations', surfID, xx);
            c = pt.Value;
        end

        function set.coverages(s, cov, norm)
            % SET.COVERAGES  Set surface coverages of the species on an interface.
            % s.coverages = (cov, norm)
            % :param s:
            %      Instance of class :mat:func:`Interface`
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
            %
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
                if sz == nsp
                    if ((m == nsp && n == 1) || (m == 1 & n == nsp))
                        callct('surf_setCoverages', surfID, cov, norm_flag);
                    else error('wrong size for coverage array');
                    end
                else
                    error('wrong size for coverage array');
                end
            elseif isa(cov, 'char')
                callct('surf_setCoveragesByName', surfID, cov);
            end
        end

        function set.siteDensity(s, d)
            % SET.SITEDENSITY  Set the site density of a phase on an interface.
            % s.siteDensity = d
            % :param s:
            %      Instance of class :mat:func:`Interface`
            % :parameter d
            %    Double site density. Unit: kmol/m^2 for surface phases,
            %    kmol/m for edge phases.

            surfID = s.tpID;
            callct('surf_setSiteDensity', surfID, d);
        end

    end
end


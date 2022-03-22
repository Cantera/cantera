classdef Interface < handle & ThermoPhase & Kinetics

    properties
        coverages
    end

    methods
        %% Interface class constructor

        function s = Interface(src, id, p1, p2, p3, p4)
            % :parameter src:
            %    CTI or CTML file containing the interface or edge phase.
            % :parameter id:
            %    Name of the interface or edge phase in the source file.
            % :parameter p1/p2/p3/p4:
            %    Adjoining phase to the interface;
            % :return:
            %    Instance of class 'Interface'.

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
            % Get the surface coverages of the species on an interface.
            %
            % :return:
            %    Vector of length "n_surf_species" for coverage.

            surfID = s.tpID;
            nsp = s.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'surf_getCoverages', surfID, xx);
            c = pt.Value;
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
            calllib(ct, 'surf_getConcentrations', surfID, xx);
            c = pt.Value;
        end

        function set.coverages(s, cov, norm)
            % Set surface coverages of the species on an interface.
            %
            % :parameter cov:
            %    Coverage of the species. "Cov" can be either a vector of
            %    length "n_surf_species", or a string in the format
            %    "Species:Coverage, Species:Coverage".

            if nargin == 3 && strcmp(norm, 'nonorm')
                norm_flag = 0;
            else
                norm_flag = 1;
            end

            surfID = s.tr_id;
            nsp = s.nSpecies;
            [m, n] = size(cov);

            if isa(cov, 'double')
                sz = length(cov);
                if sz == nsp
                    if ((m == nsp && n == 1) || (m == 1 & n == nsp))
                        calllib(ct, 'surf_setCoverages', surfID, cov, norm_flag);
                    else error('wrong size for coverage array');
                    end
                else
                    error('wrong size for coverage array');
                end
            elseif isa(cov, 'char')
                calllib(ct, 'surf_setCoveragesByName', surfID, cov);
            end
        end

    end
end

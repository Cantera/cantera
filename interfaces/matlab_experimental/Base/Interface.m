classdef Interface < handle

    properties
        th
        kin
    end

    properties(Constant = true)
        lib = 'cantera_shared'
    end

    methods
        %% Interface class constructor

        function s = Interface(src, id, p1, p2, p3, p4)
            % :param src:
            %    CTI or CTML file containing the interface or edge phase.
            % :param id:
            %    Name of the interface or edge phase in the source file.
            % :param p1/P2/P3/P4:
            %    Adjoining phase to the interface;
            % :return:
            %    Instance of class 'Interface'.

            checklib;
            t = ThermoPhase(src, id);
            if nargin == 2
                k = Kinetics(t, src, id);
            elseif nargin == 3
                k = Kinetics(t, src, id, p1);
            elseif nargin == 4
                k = Kinetics(t, src, id, p1, p2);
            elseif nargin == 5
                k = Kinetics(t, src, id, p1, p2, p3);
            elseif nargin == 6
                k = Kinetics(t, src, id, p1, p2, p3, p4);
            end

            s.kin = k;
            s.th = t;
        end

        %% Interface methods

        function c = coverages(s)
            % Get the surface coverages of the species on an interface.
            %
            % :return:
            %    If no output value is assigned, a bar graph will be
            %    plotted. Otherwise, a vector of length "n_surf_species"
            %    will be returned.

            checklib;
            surf_id = s.th.tr_id;
            nsp = s.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(s.lib, 'surf_getCoverages', surf_id, xx);
            c = pt.Value;

            if nargout == 0
                figure
                set(gcf, 'Name', 'Coverages')
                bar(c);
                colormap(summer);
                nm = speciesNames(s);
                set(gca, 'XTickLabel', nm);
                xlabel('Species Name');
                ylabel('Coverage');
                title('Surface Species Coverages');
            end
        end

        function c = concentrations(s)
            % Get the concentrations of the species on an interface.
            %
            % :return:
            %    If no output value is assigned, a bar graph will be
            %    plotted. Otherwise, a vector of length "n_surf_species"
            %    will be returned.

            checklib;
            surf_id = s.th.tr_id;
            nsp = s.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(s.lib, 'surf_getConcentrations', surf_id, xx);
            c = pt.Value;

            if nargout == 0
                figure
                set(gcf, 'Name', 'Concentrations')
                bar(c);
                colormap(summer);
                nm = speciesNames(s);
                set(gca, 'XTickLabel', nm);
                xlabel('Species Name');
                ylabel('Concentration [kmol/m^2]');
                title('Surface Species Concentrations');
            end
        end

        function setCoverages(s, cov, norm)
            % Set surface coverages of the species on an interface.
            %
            % :param cov:
            %    Coverage of the species. "Cov" can be either a vector of
            %    length "n_surf_species", or a string in the format
            %    "Species:Coverage, Species:Coverage".

            checklib;

            if nargin == 3 && strcmp(norm, 'nonorm')
                norm_flag = 0;
            else
                norm_flag = 1;
            end

            surf_id = s.th.tr_id;
            nsp = s.th.nSpecies;
            [m, n] = size(cov);

            if isa(cov, 'double')
                sz = length(cov);
                if sz == nsp
                    if ((m == nsp && n == 1) || (m == 1 & n == nsp))
                        calllib(s.lib, 'surf_setCoverages', surf_id, cov, norm_flag);
                    else error('wrong size for coverage array');
                    end
                else
                    error('wrong size for coverage array');
                end
            elseif isa(cov, 'char')
                calllib(s.lib, 'surf_setCoveragesByName', surf_id, cov);
            end
        end

    end
end

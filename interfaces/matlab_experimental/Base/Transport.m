classdef Transport < handle

    properties
        th
        tr_id
    end

    methods
        %% Transport class constructor

        function tr = Transport(tp, model, loglevel)
            checklib;
            tr.tr_id = 0;
            if nargin == 2
                model = 'default'
            end

            if nargin < 3
                loglevel = 4;
            end

            if ~isa(tp, 'ThermoPhase')
                error(['The first argument must be an ', ...
                      'instance of class ThermoPhase']);
            else
                tr.th = tp;
                if strcmp(model, 'default')
                    tr.tr_id = calllib(ct, 'trans_newDefault', ...
                                       tp.tp_id, loglevel);
                else
                    tr.tr_id = calllib(ct, 'trans_new', model, ...
                                       tp.tp_id, loglevel);
                end
            end
            tr.tp_id = tp.tp_id;
        end

        %% Utility methods

        function tr_clear(tr)
            % Delete the kernel object.

            checklib;
            calllib(ct, 'trans_del', tr.tr_id);
        end

        %% Transport Methods

        function v = viscosity(tr)
            % Get the dynamic viscosity.
            %
            % return:
            %    Double dynamic viscosity. Unit: Pa*s.

            checklib;
            v = calllib(ct, 'trans_viscosity', tr.tr_id);
            if v == -1.0
                error(geterr);
            elseif v < 0.0
                error('exception raised');
            end
        end

        function v = thermalConductivity(tr)
            % Get the thermal conductivity.
            %
            % return:
            %    Double thermal conductivity. Unit: W/m-K.

            checklib;
            v = calllib(ct, 'trans_thermalConductivity', tr.tr_id);
            if v == -1.0
                error(geterr);
            elseif v < 0.0
                error('exception raised');
            end
        end

        function v = electricalConductivity(tr)
            % Get the electrical conductivity.
            %
            % return:
            %    Double electrical conductivity. Unit: S/m.

            checklib;
            v = calllib(ct, 'trans_electricalConductivity', tr.tr_id);
            if v == -1.0
                error(geterr);
            elseif v < 0.0
                error('exception raised');
            end
        end

        function v = mixDiffCoeffs(tr)
            % Get the mixture-averaged diffusion coefficients.
            %
            % return:
            %    Vector of length nSpecies with the mixture-averaged
            %    diffusion coefficients. Unit: m^2/s.

            checklib;
            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'trans_getMixDiffCoeffs', tr.tr_id, nsp, pt);
            v = pt.Value;
        end

        function v = thermalDiffCoeffs(tr)
            % Get the thermal diffusion coefficients.
            %
            % return:
            %    Vector of length nSpecies with the thermal diffusion
            %    coefficients.

            checklib;
            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'trans_getThermalDiffCoeffs', tr.tr_id, nsp, pt);
            v = pt.Value;
        end

        function v = binDiffCoeffs(tr)
            % Get the binary diffusion coefficients.
            %
            % return:
            %    A matrix of binary diffusion coefficients. The matrix is
            %    symmetric: d(i, j) = d(j, i). Unit: m^2/s.

            checklib;
            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'trans_getBinDiffCoeffs', tr.tr_id, nsp, pt);
            v = pt.Value;
        end

        function v = multiDiffCoeffs(tr)
            % Get the multicomponent diffusion coefficients.
            %
            % return:
            %    Vector of length nSpecies with the multicomponent
            %    diffusion coefficients. Unit: m^2/s.

            checklib;
            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            calllib(ct, 'trans_getMultiDiffCoeffs', tr.tr_id, nsp, pt);
            v = pt.Value;
        end

        function setParameters(tr, type, k, p)
            % Set the parameters.
            %
            % parameter type:
            % parameter k:
            % parameter p:

            checklib;
            calllib(ct, 'trans_setParameters', tr.tr_id, type, k, p);
        end

        function setThermalConductivity(tr, lam)
            % Set the thermal conductivity.
            %
            % parameter lam:
            %    Thermal conductivity in W/(m-K).

            checklib;
            tr.setParameters(1, 0, lam);
        end

    end
end

classdef Transport < handle

    properties
        th
        trID
    end

    methods
        %% Transport class constructor

        function tr = Transport(tp, model, loglevel)
            checklib;
            tr.trID = 0;
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
                    tr.trID = callct('trans_newDefault', ...
                                       tp.tpID, loglevel);
                else
                    tr.trID = callct('trans_new', model, ...
                                       tp.tpID, loglevel);
                end
            end
            tr.tpID = tp.tpID;
        end

        %% Utility methods

        function trClear(tr)
            % Delete the kernel object.

            callct('trans_del', tr.trID);
        end

        %% Transport Methods

        function v = viscosity(tr)
            % Get the dynamic viscosity.
            %
            % :return:
            %    Double dynamic viscosity. Unit: Pa*s.

            v = callct('trans_viscosity', tr.trID);
            if v == -1.0
                error(geterr);
            elseif v < 0.0
                error('exception raised');
            end
        end

        function v = thermalConductivity(tr)
            % Get the thermal conductivity.
            %
            % :return:
            %    Double thermal conductivity. Unit: W/m-K.

            v = callct('trans_thermalConductivity', tr.trID);
            if v == -1.0
                error(geterr);
            elseif v < 0.0
                error('exception raised');
            end
        end

        function v = electricalConductivity(tr)
            % Get the electrical conductivity.
            %
            % :return:
            %    Double electrical conductivity. Unit: S/m.

            v = callct('trans_electricalConductivity', tr.trID);
            if v == -1.0
                error(geterr);
            elseif v < 0.0
                error('exception raised');
            end
        end

        function v = mixDiffCoeffs(tr)
            % Get the mixture-averaged diffusion coefficients.
            %
            % :return:
            %    Vector of length nSpecies with the mixture-averaged
            %    diffusion coefficients. Unit: m^2/s.

            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getMixDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = thermalDiffCoeffs(tr)
            % Get the thermal diffusion coefficients.
            %
            % :return:
            %    Vector of length nSpecies with the thermal diffusion
            %    coefficients.

            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getThermalDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = binDiffCoeffs(tr)
            % Get the binary diffusion coefficients.
            %
            % :return:
            %    A matrix of binary diffusion coefficients. The matrix is
            %    symmetric: d(i, j) = d(j, i). Unit: m^2/s.

            nsp = tr.th.nSpecies;
            xx = zeros(nsp, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getBinDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = multiDiffCoeffs(tr)
            % Get the multicomponent diffusion coefficients.
            %
            % :return:
            %    Vector of length nSpecies with the multicomponent
            %    diffusion coefficients. Unit: m^2/s.

            nsp = tr.th.nSpecies;
            xx = zeros(nsp, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getMultiDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function setParameters(tr, type, k, p)
            % Set the parameters.
            %
            % :parameter type:
            % :parameter k:
            % :parameter p:

            callct('trans_setParameters', tr.trID, type, k, p);
        end

        function setThermalConductivity(tr, lam)
            % Set the thermal conductivity.
            %
            % :parameter lam:
            %    Thermal conductivity in W/(m-K).

            tr.setParameters(1, 0, lam);
        end

    end
end

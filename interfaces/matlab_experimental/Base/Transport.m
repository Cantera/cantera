classdef Transport < handle

    properties (Access = private)
        th
    end

    properties (SetAccess = immutable)
        trID % ID of Transport object
    end

    properties (SetAccess = protected)

        viscosity % Dynamic viscosity. Unit: Pa*s.

        thermalConductivity % Thermal conductivity. Unit: W/m-K.

        electricalConductivity % Electrical conductivity. Unit: S/m.

        mixDiffCoeffs % Mixture-averaged diffusion coefficients. Unit: m^2/s.

        thermalDiffCoeffs % Thermal diffusion coefficients.

        binDiffCoeffs % Binary diffusion coefficients. Unit: m^2/s.

        multiDiffCoeffs % Multicomponent diffusion coefficients. Unit: m^2/s.

    end

    methods
        %% Transport Class Constructor

        function tr = Transport(tp, model, loglevel)
            % Transport Class
            %
            % tr = Transport(r, th, model, loglevel)
            %
            % Create a new instance of class :mat:class:`Transport`. One to three arguments
            % may be supplied. The first must be an instance of class
            % :mat:class:`ThermoPhase`. The second (optional) argument is the type of
            % model desired, specified by the string ``'default'``, ``'Mix'`` or
            % ``'Multi'``. ``'default'`` uses the default transport specified in the
            % phase definition. The third argument is the logging level desired.
            %
            % :param th:
            %     Instance of class :mat:class:`ThermoPhase`
            % :param model:
            %     String indicating the transport model to use. Possible values
            %     are ``'default'``, ``'None'``, ``'Mix'``, and ``'Multi'``.
            %     Optional.
            % :param loglevel:
            %     Level of diagnostic logging. Default if not specified is 4.
            % :return:
            %     Instance of class :mat:class:`Transport`
            %
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
            end

            tr.th = tp;

            if strcmp(model, 'default')
                tr.trID = callct('trans_newDefault', tp.tpID, loglevel);
            else
                tr.trID = callct('trans_new', model, tp.tpID, loglevel);
            end

            tr.tpID = tp.tpID;
        end

        %% Transport Class Destructor

        function delete(tr)
            % Delete the kernel object.

            callct('trans_del', tr.trID);
        end

        %% Transport Get Methods

        function v = get.viscosity(tr)
            v = callct('trans_viscosity', tr.trID);
        end

        function v = get.thermalConductivity(tr)
            v = callct('trans_thermalConductivity', tr.trID);
        end

        function v = get.electricalConductivity(tr)
            v = callct('trans_electricalConductivity', tr.trID);
        end

        function v = get.mixDiffCoeffs(tr)
            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getMixDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = get.thermalDiffCoeffs(tr)
            nsp = tr.th.nSpecies;
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getThermalDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = get.binDiffCoeffs(tr)
            nsp = tr.th.nSpecies;
            xx = zeros(nsp, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getBinDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = get.multiDiffCoeffs(tr)
            nsp = tr.th.nSpecies;
            xx = zeros(nsp, nsp);
            pt = libpointer('doublePtr', xx);
            callct('trans_getMultiDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        %% Transport Set Methods

        function setParameters(tr, type, k, p)
            % Set the parameters.
            %
            % tr.setParameters(type, k, p)
            %
            % :param type:
            % :param k:
            % :param p:

            callct('trans_setParameters', tr.trID, type, k, p);
        end

        function setThermalConductivity(tr, lam)
            % Set the thermal conductivity.
            %
            % tr.setThermalConductivity(lam)
            %
            % :param lam:
            %    Thermal conductivity in W/(m-K).
            %

            tr.setParameters(1, 0, lam);
        end

    end

end

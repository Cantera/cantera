classdef Transport < handle
    % Transport Class ::
    %
    %     >> tr = Transport(th, model, loglevel)
    %
    % Create a new instance of class :mat:class:`Transport`. One to three
    % arguments may be supplied. The first must be an instance of class
    % :mat:class:`ThermoPhase`. The second (optional) argument is the type of
    % model desired, specified by the string ``'default'``,
    % ``'mixture-averaged'`` or ``'multicomponent'``. ``'default'``
    % uses the default transport specified in the phase definition.
    % The third argument is the logging level desired.
    %
    % :param tp:
    %     Instance of class :mat:class:`ThermoPhase` OR string "clib"
    %     when called by the class constructors of :mat:class:`Solution` or
    %     :mat:class:`Interface`.
    % :param model:
    %     String indicating the transport model to use. Possible values
    %     are ``'default'``, ``'none'``, ``'mixture-averaged'``,
    %     and ``'multicomponent'``. Optional.
    % :param loglevel:
    %     Level of diagnostic logging. Default if not specified is 4. Optional.
    % :return:
    %     Instance of class :mat:class:`Transport`.

    properties (Access = protected)
        th % ID of the ThermoPhase object used to create the Transport object.
    end

    properties (SetAccess = immutable)
        trID % ID of Transport object.
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

        function tr = Transport(varargin)

            ctIsLoaded;
            tr.trID = 0;

            tp = varargin{1};

            if ischar(tp) & isnumeric(varargin{2})
                if strcmp(tp, 'clib')
                    tr.trID = ctFunc('soln_transport', varargin{2});
                    return
                end
            end

            if nargin < 2
                model = 'default'
            else
                model = varargin{2};
            end

            if nargin < 3
                loglevel = 4;
            end

            if ~isa(tp, 'ThermoPhase')
                error(['The first argument must be an ', ...
                       'instance of class ThermoPhase']);
            end

            tr.th = tp.tpID;

            if strcmp(model, 'default')
                tr.trID = ctFunc('trans_newDefault', tp.tpID, loglevel);
            else
                tr.trID = ctFunc('trans_new', model, tp.tpID, loglevel);
            end

        end

        %% Transport Destructor Methods

        function delete(tr)
            % Delete the :mat:class:`Transport` object.

            if ~isa(tr, 'Solution')
                ctFunc('trans_del', tr.trID);
            end
        end

        %% Transport Get Methods

        function v = get.viscosity(tr)
            v = ctFunc('trans_viscosity', tr.trID);
        end

        function v = get.thermalConductivity(tr)
            v = ctFunc('trans_thermalConductivity', tr.trID);
        end

        function v = get.electricalConductivity(tr)
            v = ctFunc('trans_electricalConductivity', tr.trID);
        end

        function v = get.mixDiffCoeffs(tr)
            nsp = ctFunc('thermo_nSpecies', tr.th);
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('trans_getMixDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = get.thermalDiffCoeffs(tr)
            nsp = ctFunc('thermo_nSpecies', tr.th);
            xx = zeros(1, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('trans_getThermalDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = get.binDiffCoeffs(tr)
            nsp = ctFunc('thermo_nSpecies', tr.th);
            xx = zeros(nsp, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('trans_getBinDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

        function v = get.multiDiffCoeffs(tr)
            nsp = ctFunc('thermo_nSpecies', tr.th);
            xx = zeros(nsp, nsp);
            pt = libpointer('doublePtr', xx);
            ctFunc('trans_getMultiDiffCoeffs', tr.trID, nsp, pt);
            v = pt.Value;
        end

    end

end

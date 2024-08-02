classdef Transport < handle
    % Transport Class ::
    %
    %     >> tr = Transport(id)
    %
    % Retrieve instance of class :mat:class:`Transport` associated with a
    % :mat:class:`Solution` object.
    %
    % :param id:
    %     Integer ID of the solution holding the :mat:class:`Transport` object.
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

        function tr = Transport(id)
            % Create a :mat:class:`Transport` object.
            ctIsLoaded;

            if ~isnumeric(id)
                error('Invalid argument: constructor requires integer solution ID.')
            end

            tr.trID = ctFunc('soln_transport', id);
            tr.th = ctFunc('soln_thermo', id);
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

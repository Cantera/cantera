classdef Transport < handle
    % Transport Class ::
    %
    %     >> tr = Transport(id)
    %
    % Retrieve instance of class :mat:class:`Transport` associated with a
    % :mat:class:`Solution` object. The constructor is called whenever a new
    % :mat:class:`Solution` is instantiated and should not be used directly.
    %
    % :param id:
    %     Integer ID of the solution holding the :mat:class:`Transport` object.
    % :return:
    %     Instance of class :mat:class:`Transport`.

    properties (Access = protected)
        th % ID of the ThermoPhase object used to create the Transport object.
    end

    properties (SetAccess = protected)

        trID % ID of :mat:class:`Transport` object.

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

            tr.trID = ctFunc('sol_transport', id);
            tr.th = ctFunc('sol_thermo', id);
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
            v = ctArray('trans_getMixDiffCoeffs', nsp, tr.trID);
        end

        function v = get.thermalDiffCoeffs(tr)
            nsp = ctFunc('thermo_nSpecies', tr.th);
            v = ctArray('trans_getThermalDiffCoeffs', nsp, tr.trID);
        end

        function v = get.binDiffCoeffs(tr)
            nsp = ctFunc('thermo_nSpecies', tr.th);
            val = ctArray('trans_getBinaryDiffCoeffs', nsp*nsp, tr.trID, nsp);
            v = reshape(val, [nsp, nsp]);
        end

        function v = get.multiDiffCoeffs(tr)
            nsp = ctFunc('thermo_nSpecies', tr.th);
            val = ctArray('trans_getMultiDiffCoeffs', nsp*nsp, tr.trID, nsp);
            v = reshape(val, [nsp, nsp]);
        end

    end

end

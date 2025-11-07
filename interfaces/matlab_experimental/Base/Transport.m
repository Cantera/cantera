classdef (Abstract) Transport < handle
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

    properties (Access = protected)
        th % ID of the ThermoPhase object used to create the Transport object.
    end

    properties (SetAccess = protected)

        trID % ID of :mat:class:`Transport` object.

        viscosity % Dynamic viscosity [Pa·s].

        thermalConductivity % Thermal conductivity [W/m/K]

        electricalConductivity % Electrical conductivity [S/m].

        mixDiffCoeffs % Mixture-averaged diffusion coefficients. [m²/s].

        thermalDiffCoeffs % Thermal diffusion coefficients [kg/m/s].

        binDiffCoeffs % Binary diffusion coefficients [m²/s].

        multiDiffCoeffs % Multicomponent diffusion coefficients [m²/s].

    end

    methods
        %% Transport Class Constructor

        function tr = Transport(id)
            arguments
                id (1,1) double {mustBeInteger}
            end

            tr.trID = ctFunc('mSol_transport', id);
            tr.th = ctFunc('mSol_thermo', id);
        end

        %% Transport Get Methods

        function v = get.viscosity(obj)
            v = ctFunc('mTrans_viscosity', obj.trID);
        end

        function v = get.thermalConductivity(obj)
            v = ctFunc('mTrans_thermalConductivity', obj.trID);
        end

        function v = get.electricalConductivity(obj)
            v = ctFunc('mTrans_electricalConductivity', obj.trID);
        end

        function v = get.mixDiffCoeffs(obj)
            nsp = ctFunc('mThermo_nSpecies', obj.th);
            v = ctArray('mTrans_getMixDiffCoeffs', nsp, obj.trID);
        end

        function v = get.thermalDiffCoeffs(obj)
            nsp = ctFunc('mThermo_nSpecies', obj.th);
            v = ctArray('mTrans_getThermalDiffCoeffs', nsp, obj.trID);
        end

        function v = get.binDiffCoeffs(obj)
            nsp = ctFunc('mThermo_nSpecies', obj.th);
            val = ctArray('mTrans_getBinaryDiffCoeffs', nsp*nsp, obj.trID, nsp);
            v = reshape(val, [nsp, nsp]);
        end

        function v = get.multiDiffCoeffs(obj)
            nsp = ctFunc('mThermo_nSpecies', obj.th);
            val = ctArray('mTrans_getMultiDiffCoeffs', nsp*nsp, obj.trID, nsp);
            v = reshape(val, [nsp, nsp]);
        end

    end

end

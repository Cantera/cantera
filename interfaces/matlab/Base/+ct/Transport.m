classdef (Abstract) Transport < handle
    % Transport Class ::
    %
    %     >> tr = ct.Transport(id)
    %
    % Retrieve instance of class :mat:class:`ct.Transport` associated with a
    % :mat:class:`ct.Solution` object. The constructor is called whenever a new
    % :mat:class:`ct.Solution` is instantiated and should not be used directly.
    %
    % :param id:
    %     Integer ID of the solution holding the :mat:class:`ct.Transport` object.

    properties (Access = protected)
        th % ID of the ThermoPhase object used to create the Transport object.
    end

    properties (SetAccess = protected)

        trID % ID of :mat:class:`ct.Transport` object.

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

        function obj = Transport(id)
            arguments
                id (1,1) double {mustBeInteger}
            end

            obj.trID = ct.impl.call('mSol_transport', id);
            obj.th = ct.impl.call('mSol_thermo', id);
        end

        %% Transport Get Methods

        function v = get.viscosity(obj)
            v = ct.impl.call('mTrans_viscosity', obj.trID);
        end

        function v = get.thermalConductivity(obj)
            v = ct.impl.call('mTrans_thermalConductivity', obj.trID);
        end

        function v = get.electricalConductivity(obj)
            v = ct.impl.call('mTrans_electricalConductivity', obj.trID);
        end

        function v = get.mixDiffCoeffs(obj)
            nsp = ct.impl.call('mThermo_nSpecies', obj.th);
            v = ct.impl.getArray('mTrans_getMixDiffCoeffs', nsp, obj.trID);
        end

        function v = get.thermalDiffCoeffs(obj)
            nsp = ct.impl.call('mThermo_nSpecies', obj.th);
            v = ct.impl.getArray('mTrans_getThermalDiffCoeffs', nsp, obj.trID);
        end

        function v = get.binDiffCoeffs(obj)
            nsp = ct.impl.call('mThermo_nSpecies', obj.th);
            val = ct.impl.getArray('mTrans_getBinaryDiffCoeffs', nsp*nsp, obj.trID, nsp);
            v = reshape(val, [nsp, nsp]);
        end

        function v = get.multiDiffCoeffs(obj)
            nsp = ct.impl.call('mThermo_nSpecies', obj.th);
            val = ct.impl.getArray('mTrans_getMultiDiffCoeffs', nsp*nsp, obj.trID, nsp);
            v = reshape(val, [nsp, nsp]);
        end

    end

end

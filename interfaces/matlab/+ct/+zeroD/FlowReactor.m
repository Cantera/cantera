classdef FlowReactor < ct.zeroD.ReactorBase
    % Create a flow reactor object. ::
    %
    %     >> r = ct.zeroD.FlowReactor(phase, name, clone)
    %
    % A reactor representing adiabatic plug flow in a constant-area
    % duct. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = ct.zeroD.FlowReactor(gas)    % a reactor containing a gas
    %
    % See also: :mat:class:`ct.zeroD.ReactorBase`
    %
    % :param phase:
    %     Cantera :mat:class:`ct.Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `phase` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    %    (optional; default is true).

    properties (SetAccess = public)

        % Mass flow rate through the reactor [kg/s].
        %
        % Setting the mass flow rate sets the flow speed based on the current
        % density of the reactor contents and the reactor area. The value read
        % back is computed from the flow speed, density and area at the end of
        % the last call to "advance" or "step".
        massFlowRate

    end

    methods

        function obj = FlowReactor(phase, name, clone)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            id = ct.impl.call('mReactor_new', 'FlowReactor', phase.solnID, clone, name);
            obj@ct.zeroD.ReactorBase(id);
        end

        %% FlowReactor Get Methods

        function rate = get.massFlowRate(obj)
            rate = ct.impl.call('mReactor_massFlowRate', obj.id);
        end

        %% FlowReactor Set Methods

        function set.massFlowRate(obj, MFR)
            ct.impl.call('mReactor_setMassFlowRate', obj.id, MFR);
        end

    end
end

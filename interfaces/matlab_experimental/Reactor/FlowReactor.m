classdef FlowReactor < ReactorBase
    % Create a flow reactor object. ::
    %
    %     >> r = FlowReactor(phase, name, clone)
    %
    % A reactor representing adiabatic plug flow in a constant-area
    % duct. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = FlowReactor(gas)    % a reactor containing a gas
    %
    % See also: :mat:class:`ReactorBase`
    %
    % :param phase:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `phase` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    %    (optional; default is true).

    properties (SetAccess = public)

        massFlowRate % Mass flow rate in kg/s.

    end

    methods

        function r = FlowReactor(phase, name, clone)
            arguments
                phase (1,1) Solution
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            id = ctFunc('mReactor_new', 'FlowReactor', phase.solnID, clone, name);
            r@ReactorBase(id);
        end

        %% FlowReactor Get Methods

        function flag = get.massFlowRate(r)
            rate = ctFunc('mReactor_massFlowRate', r.id);
        end

        %% FlowReactor Set Methods

        function set.massFlowRate(r, MFR)
            ctFunc('mReactor_setMassFlowRate', r.id, MFR);
            r.massFlowRate = MFR;
        end

    end
end

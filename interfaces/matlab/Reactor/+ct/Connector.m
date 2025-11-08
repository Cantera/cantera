classdef (Abstract) Connector < handle
    % Connector Class ::
    %
    %     >> c = Connector(typ, r1, r2, name)
    %
    % Base class for walls and flow devices.
    %
    % See also: :mat:class:`FlowDevice`, :mat:class:`Wall`
    %
    % :param typ:
    %     Type of connector.
    % :param r1:
    %     Reactor one.
    % :param r2:
    %     Reactor two.
    % :param name:
    %     Connector name (optional; default is ``(none)``).

    properties (SetAccess = immutable)

        id = -1  % ID of Connector object.

    end

    properties (SetAccess = public)

        type  % Name of connector.

        name  % Name of connector.

    end

    methods
        %% Connector Class Constructor

        function c = Connector(typ, r1, r2, name)
            arguments
                typ (1,1) string
                r1 (1,1) ct.ReactorBase
                r2 (1,1) ct.ReactorBase
                name (1,1) string = "(none)"
            end

            c.id = ct.impl.call('mConnector_new', typ, r1.id, r2.id, name);
        end

        %% Connector Class Destructor

        function delete(obj)
            % Delete the :mat:class:`Connector` object.
            if obj.id >= 0
                ct.impl.call('mConnector_del', obj.id);
            end
        end

        %% Connector Get Methods

        function typ = get.type(obj)
            typ = ct.impl.getString('mConnector_type', obj.id);
        end

        function name = get.name(obj)
            name = ct.impl.getString('mConnector_name', obj.id);
        end

        %% Connector Set Methods

        function set.name(obj, name)
            ct.impl.call('mConnector_setName', obj.id, name);
        end

    end

end

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
    % :return:
    %     Instance of class :mat:class:`Connector`.

    properties (SetAccess = immutable)

        id % ID of Connector object.

    end

    properties (SetAccess = public)

        type  % Name of connector.

        name  % Name of connector.

    end

    methods
        %% Connector Class Constructor

        function c = Connector(typ, r1, r2, name)
            % Create a :mat:class:`Connector` object.

            arguments
                typ (1,1) string
                r1 {mustBeA(r1, 'ReactorBase')}
                r2 {mustBeA(r2, 'ReactorBase')}
                name (1,1) string = "(none)"
            end

            ctIsLoaded;

            c.id = ctFunc('mConnector_new', typ, r1.id, r2.id, name);
        end

        %% Connector Class Destructor

        function delete(c)
            % Delete the :mat:class:`Connector` object.
            ctFunc('mConnector_del', c.id);
        end

        %% Connector Get Methods

        function typ = get.type(c)
            typ = ctString('mConnector_type', c.id);
        end

        function name = get.name(c)
            name = ctString('mConnector_name', c.id);
        end

        %% Connector Set Methods

        function set.name(c, name)
            ctFunc('mConnector_setName', c.id, name);
        end

    end

end

classdef Connector < handle
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

            ctIsLoaded;

            if nargin < 3
                error('please specify type and reactors');
            end
            if nargin < 4
                name = '(none)';
            end

            if ~isa(r1, 'Reactor') || ~isa(r1, 'Reactor')
                error(['Connectors can only be installed between', ...
                       'reactors or reservoirs']);
            end

            c.id = ctFunc('connector_new', typ, r1.id, r2.id, name);
        end

        %% Connector Class Destructor

        function delete(c)
            % Delete the :mat:class:`Connector` object.
            if isempty(c.id)
                return
            end
            ctFunc('connector_del', c.id);
        end

        %% Connector Get Methods

        function typ = get.type(c)
            typ = ctString('connector_type', c.id);
        end

        function name = get.name(c)
            name = ctString('connector_name', c.id);
        end

        %% Connector Set Methods

        function set.name(c, name)
            ctFunc('connector_setName', c.id, name);
        end

    end

end

classdef Reservoir < Reactor
    % Create a :mat:class:`Reservoir` object. ::
    %
    %     >> r = Reservoir(contents, name)
    %
    % A :mat:class:`Reservoir` is an instance of class :mat:class:`Reactor`
    % configured so that its intensive state is constant in time. A reservoir
    % may be thought of as infinite in extent, perfectly mixed,
    % and non-reacting, so that fluid may be extracted or added without
    % changing the composition or thermodynamic state. Note that even
    % if the reaction mechanism associated with the fluid in the
    % reactor defines reactions, they are disabled within
    % reservoirs. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = Reservoir(gas)    % a reservoir containing a gas
    %
    % See also: :mat:class:`Reactor`
    %
    % :param contents:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reservoir name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`Reactor`.

    methods

        function r = Reservoir(contents, name)
            % Constructor

            if nargin < 2
                name = '(none)';
            end

            r@Reactor(contents, 'Reservoir', name);
        end

    end
end

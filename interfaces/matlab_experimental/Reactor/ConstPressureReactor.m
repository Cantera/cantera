classdef ConstPressureReactor < Reactor
    % Create a constant pressure reactor object. ::
    %
    %     >> r = ConstPressureReactor(contents, name)
    %
    % A :mat:class:`ConstPressureReactor` is an instance of class
    % :mat:class:`Reactor` where the pressure is held constant. The volume
    % is not a state variable, but instead takes on whatever value is
    % consistent with holding the pressure constant. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = ConstPressureReactor(contents)    % a reactor containing contents
    %
    % See also: :mat:class:`Reactor`
    %
    % :param contents:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`ConstPressureReactor`.

    methods

        function r = ConstPressureReactor(contents, name)
            % Constructor

            if nargin < 2
                name = '(none)';
            end

            r@Reactor(contents, 'ConstPressureReactor', name);
        end

    end
end

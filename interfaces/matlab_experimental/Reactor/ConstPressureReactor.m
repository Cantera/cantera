classdef ConstPressureReactor < Reactor
    % Create a constant pressure reactor object.
    %
    % r = ConstPressureReactor(contents)
    %
    % A :mat:class:`ConstPressureReactor` is an instance of class
    % :mat:class:`Reactor` where the pressure is held constant. The volume
    % is not a state variable, but instead takes on whatever value is
    % consistent with holding the pressure constant. Examples:
    %
    % .. code-block:: matlab
    %
    %     r1 = ConstPressureReactor         % an empty reactor
    %     r2 = ConstPressureReactor(contents)    % a reactor containing contents
    %
    % See also: :mat:class:`Reactor`
    %
    % :param contents:
    %     Cantera :mat:class:`Solution` to be set as the contents of the
    %     reactor
    % :return:
    %     Instance of class :mat:class:`ConstPressureReactor`

    methods

        % Constructor
        function r = ConstPressureReactor(contents)
            if nargin == 0
                contents = 0;
            end

            r = r@Reactor(contents, 'ConstPressureReactor');
        end

    end
end

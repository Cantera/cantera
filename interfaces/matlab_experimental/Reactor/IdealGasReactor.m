classdef IdealGasReactor < Reactor
    % Create a reactor with an ideal gas. ::
    %
    %     >> r = IdealGasReactor(contents)
    %
    % An IdealGasReactor is an instance of class Reactor where the governing
    % equations are specialized for the ideal gas equation of state (and do not
    % work correctly with other thermodynamic models). Examples:
    %
    % .. code-block:: matlab
    %
    %     r1 = IdealGasReactor         % an empty reactor
    %     r2 = IdealGasReactor(gas)    % a reactor containing a gas
    %
    % See also: :mat:class:`Reactor`
    %
    % :param contents:
    %     Cantera :mat:class:`Solution` to be set as the contents of the
    %     reactor
    % :return:
    %     Instance of class :mat:class:`IdealGasReactor`

    methods

        function r = IdealGasReactor(contents)
            % Constructor

            if nargin == 0
                contents = 0;
            end

            r@Reactor(contents, 'IdealGasReactor');
        end

    end

end

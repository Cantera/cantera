classdef fplus < Func
    % Fplus - Get a functor representing the sum of two input functors. ::
    %
    %     >> f = fplus(a, b)
    %
    % :param a:
    %     Instance of class :mat:class:`Func`
    % :param b:
    %     Instance of class :mat:class:`Func`
    % :return:
    %     Instance of class :mat:class:`fplus`
    %

    methods

        function f = fplus(a, b)
            % Constructor

            f@Func('sum', a, b);
        end

    end

end

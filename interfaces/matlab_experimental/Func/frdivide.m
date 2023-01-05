classdef frdivide < Func
    % Frdivide - Get a functor representing the ratio of two input functors. ::
    %
    %     >> f = frdivide(a, b)
    %
    % :param a:
    %     Instance of class :mat:class:`Func`
    % :param b:
    %     Instance of class :mat:class:`Func`
    % :return:
    %     Instance of class :mat:class:`frdivide`
    %
    methods

        function f = frdivide(a, b)
            % Constructor

            f@Func('ratio', a, b);
        end

    end

end

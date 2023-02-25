classdef ftimes < Func
    % Get a functor representing the product of two input functors. ::
    %
    %     >> f = ftimes(a, b)
    %
    % :param a:
    %     Instance of class :mat:class:`Func`.
    % :param b:
    %     Instance of class :mat:class:`Func`.
    % :return:
    %     Instance of class :mat:class:`ftimes`.

    methods

        function f = ftimes(a, b)
            % Constructor

            f@Func('prod', a, b);
        end

    end

end

function r = fplus(a, b)
    % Get a functor representing the sum of two input functors.
    %
    % r = fplus(a, b)
    %
    % :param a:
    %     Instance of class :mat:func:`Func`
    % :param b:
    %     Instance of class :mat:func:`Func`
    % :return:
    %     Instance of class :mat:func:`Func`
    %
    r = Func('sum', a, b);
end

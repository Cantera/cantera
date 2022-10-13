function r = ftimes(a, b)
     % Create a functor that multiplies two other functors.
     %
     % r = ftimes(a, b)
     %
     % :param a:
     %     Instance of class :mat:func:`Func`
     % :param b:
     %     Instance of class :mat:func:`Func`
     % :return:
     %     Instance of class :mat:func:`Func`
     %
     r = Func('prod', a, b);
end

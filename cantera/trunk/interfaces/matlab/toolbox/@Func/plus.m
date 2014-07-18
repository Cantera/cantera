function r = plus(a, b)
% PLUS  Get a functor representing the sum of two input functors.
% r = plus(a, b)
% :param a:
%     Instance of class :mat:func:`Func`
% :param b:
%     Instance of class :mat:func:`Func`
% :return:
%     Instance of class :mat:func:`Func`
%
% PLUS - Return a functor representing the sum of two functors a
% and b.
%
r = Func('sum', a, b);

function r = rdivide(a, b)
% RDIVIDE  Get a functor that is the ratio of the input functors.
% r = rdivide(a,b)
% :param a:
%     Instance of class :mat:func:`Func`
% :param b:
%     Instance of class :mat:func:`Func`
% :return:
%     Instance of class :mat:func:`Func`
%

r = Func('ratio', a, b);

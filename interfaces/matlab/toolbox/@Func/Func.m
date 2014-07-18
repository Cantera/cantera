function x = Func(typ, n, p)
% FUNC  Func class constructor.
% x = Func(typ, n, p)
% A class for functors.
% A functor is an object that behaves like a function. Cantera
% defines a set of functors to use to create arbitrary functions to
% specify things like heat fluxes, piston speeds, etc., in reactor
% network simulations. Of course, they can be used for other things
% too.
%
% The main feature of a functor class is that it overloads the ``()``
% operator to evaluate the function. For example, suppose object
% ``f`` is a functor that evaluates the polynomial :math:`2x^2 - 3x + 1`.
% Then writing ``f(2)`` would cause the method that evaluates the
% function to be invoked, and would pass it the argument ``2``. The
% return value would of course be 3.
%
% The types of functors you can create in Cantera are these:
%
% 1. A polynomial
% 2. A Fourier series
% 3. A sum of Arrhenius terms
% 4. A Gaussian.
%
% You can also create composite functors by adding, multiplying, or
% dividing these basic functors, or other composite functors.
%
% Note: this MATLAB class shadows the underlying C++ Cantera class
% "Func1". See the Cantera C++ documentation for more details.
%
% See also: :mat:func:`polynom`, :mat:func:`gaussian`, :mat:func:`plus`,
% :mat:func:`rdivide`, :mat:func:`times`
%
% :param typ:
%     String indicating type of functor to create. Possible values are:
%
%     * ``'polynomial'``
%     * ``'fourier'``
%     * ``'gaussian'``
%     * ``'arrhenius'``
%     * ``'sum'``
%     * ``'diff'``
%     * ``'ratio'``
%     * ``'composite'``
%     * ``'periodic'``
%
% :param n:
%     Number of parameters required for the functor
% :param p:
%     Vector of parameters
% :return:
%     Instance of class :mat:func:`Func`

if ~isa(typ, 'char')
    error('Function type must be a string')
end

x.f1 = 0;
x.f2 = 0;
x.coeffs = 0;

itype = -1;
if strcmp(typ, 'polynomial')
    itype = 2;
elseif strcmp(typ, 'fourier')
    itype = 1;
elseif strcmp(typ, 'arrhenius')
    itype = 3;
elseif strcmp(typ, 'gaussian')
    itype = 4;
end

if itype > 0
    x.coeffs = p;
    x.index = funcmethods(0, itype, n, p);
elseif strcmp(typ, 'periodic')
    itype = 50;
    x.f1 = n;
    x.coeffs = p;
    x.index = funcmethods(0, itype, n.index, p);
else
    if strcmp(typ, 'sum')
        itype = 20;
    elseif strcmp(typ, 'diff')
        itype = 25;
    elseif strcmp(typ, 'prod')
        itype = 30;
    elseif strcmp(typ, 'ratio')
        itype = 40;
    elseif strcmp(typ, 'composite')
        itype = 60;
    end
    x.f1 = n;
    x.f2 = p;
    x.index = funcmethods(0, itype, n.index, p.index);
end

x.typ = typ;
x = class(x, 'Func');

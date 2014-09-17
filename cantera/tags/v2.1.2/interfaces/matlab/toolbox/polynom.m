function poly = polynom(coeffs)
%
% POLY - create an instance of class 'Func' representing a polynomial.
%
%   The polynomial coefficients are specified by a one-dimensional
%   array [a0 a1 .... aN].
%
%       polynom([-2 6 3])           3x^2 + 6x - 2
%       polynom([1.0 -2.5 0 0 2])   2x^4 - 2.5x + 1
%
[n m] = size(coeffs);
if n == 1
    poly = Func('polynomial',m - 1,coeffs);
elseif m == 1
    poly = Func('polynomial',n - 1,coeffs);
else
    error('wrong shape for coefficient array');
end

function poly = polynom(coeffs)
% POLY - create a polynomial Func instance
%   
[n m] = size(coeffs);
if n == 1
  poly = Func('polynomial',m - 1,coeffs)
elseif m == 1
  poly = Func('polynomial',n - 1,coeffs)
else
  error('wrong shape for coefficient array')
end


function g = gaussian(A, x0, FWHM)
% POLY - create a Gaussian Func instance
%        
%       gaussian(A, x0, FWHM)
%       
%          A     - value at x = x0
%          x0    - location of maximum
%          FWHM  - full width at half-maximum
%
g = Func('gaussian', 0, [A x0 FWHM]);




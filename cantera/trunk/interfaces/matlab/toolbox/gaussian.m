function g = gaussian(peak, center, width)
%
% GAUSSIAN - create a Gaussian Func instance
%
%       gaussian(peak, center, width)
%
%          peak      - the peak value
%          center    - value of x at which the peak is located
%          width     - full width at half-maximum. The value of the
%                      function at center +/- (width)/2 is one-half
%                      the peak value.
%
g = Func('gaussian', 0, [peak, center, width]);

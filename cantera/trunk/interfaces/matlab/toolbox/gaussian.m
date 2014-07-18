function g = gaussian(peak, center, width)
% GAUSSIAN  Create a Gaussian :mat:func:`Func` instance.
% g = gaussian(peak, center, width)
% :param peak:
%     The peak value
% :param center:
%     Value of x at which the peak is located
% :param width:
%     Full width at half-maximum. The value of the
%     function at center +/- (width)/2 is one-half
%     the peak value.
%

g = Func('gaussian', 0, [peak, center, width]);

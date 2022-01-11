function g = gaussian(peak, center, width)
    % Create a Gaussian functor instance.
    % :param peak:
    %    The peak value.
    % :param center:
    %    Value of x at which the peak is located.
    % :param width:
    %    Full width at half-maximum. The value of the function at center
    %    +/- (width)/2 is one-half the peak value.

    g = Func('gaussian', 0, [peak, center, width]);
end

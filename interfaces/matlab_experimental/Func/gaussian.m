classdef gaussian < Func
    % Gaussian - Create a Gaussian :mat:class:`Func` instance.
    %
    % f = gaussian(peak, center, width)
    %
    % :param peak:
    %     The peak value
    % :param center:
    %     Value of x at which the peak is located
    % :param width:
    %     Full width at half-maximum. The value of the
    %     function at center +/- (width)/2 is one-half
    %     the peak value
    % :return:
    %     Instance of class :mat:class:`gaussian`
    %
    methods

        % Constructor
        function f = gaussian(peak, center, width)
            f = f@Func('gaussian', 0, [peak, center, width]);
        end

    end

end

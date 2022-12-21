classdef polynom < Func
    % Polynom - Create a polynomial :mat:class:`Func` instance.
    %
    % poly = polynom(coeffs)
    %
    % The polynomial coefficients are specified by a vector
    % ``[a0 a1 .... aN]``. Examples:
    %
    % .. code-block:: matlab
    %
    %     polynom([-2 6 3])          %3x^2 + 6.0x - 2
    %     polynom([1.0 -2.5 0 0 2])  %2x^4 - 2.5x + 1
    %
    % :param coeffs:
    %     Vector of polynomial coefficients
    % :return:
    %     Instance of class :mat:class:`polynom`
    %
    methods

        % Constructor
        function f = polynom(coeffs)
            [n m] = size(coeffs);

    [n m] = size(coeffs);

    if n == 1
        poly = Func('polynomial', m - 1, coeffs);
    elseif m == 1
        poly = Func('polynomial', n - 1, coeffs);
    else
        error('wrong shape for coefficient array');
    end

end

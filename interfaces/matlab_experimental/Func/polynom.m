classdef polynom < Func
    % Polynom - Create a polynomial :mat:class:`Func` instance. ::
    %
    %     >> f = polynom(coeffs)
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

        function f = polynom(coeffs)
            % Constructor

            [n m] = size(coeffs);

            if n == 1
                p = m - 1;
            elseif m == 1
                p = n - 1;
            else
                error('wrong shape for coefficient array');
            end

            f@Func('polynomial', p, coeffs);
        end

    end

end

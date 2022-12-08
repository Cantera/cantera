classdef Func < handle
    % Func Class
    %
    % x = Func(typ, n, p)
    %
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
    % :param n:
    %     Number of parameters required for the functor
    % :param p:
    %     Vector of parameters
    % :return:
    %     Instance of class :mat:func:`Func`
    %

    properties (SetAccess = immutable)
        f1
        f2
        coeffs
        id
        typ
    end

    methods
        %% Func Class Constructor

        function x = Func(typ, n, p)
            checklib;

            if ~isa(typ, 'char')
                error('Function type must be a string');
            end

            x.f1 = 0;
            x.f2 = 0;
            x.coeffs = 0;
            itype = -1;

            function nn = newFunc(itype, n, p)
                % helper function to pass the correct parameters to the C
                % library
                if itype < 20
                    [m, n] = size(p);
                    lenp = m * n;
                    nn = callct('func_new', itype, n, lenp, p);
                elseif itype < 45
                    m = n;
                    nn = callct('func_new', itype, n, m, 0);
                else
                    nn = callct('func_new', itype, n, 0, p);
                end

            end

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
                x.id = newFunc(itype, n, p);
            elseif strcmp(typ, 'periodic')
                itype = 50;
                x.f1 = n;
                x.coeffs = p;
                x.id = newFunc(itype, n.id, p);
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
                x.id = newFunc(itype, n.id, p.id);
            end

            x.typ = typ;
        end

        %% Func Class Destructor

        function delete(f)
            % Delete the C++ Func1 object.

            callct('func_del', f.id);
        end

        %% Func Class Utility Methods

        function display(f)
            % Display the equation of the input function on the terminal.

            disp(' ');
            disp([inputname(1), ' = '])
            disp(' ');
            disp(['   ' f.char])
            disp(' ');
        end

        function b = subsref(a, s)
            % Redefine subscripted references for functors.
            %
            % b = a.subsref(s)
            %
            % :param a:
            %     Instance of class :mat:func:`Func`
            % :param s:
            %     Value at which the function should be evaluated.
            % :return:
            %     Returns the value of the function evaluated at ``s``
            %
            if strcmp(s.type, '()')
                ind = s.subs{:};
                b = zeros(1, length(ind));

                for k = 1:length(ind)
                    b(k) = callct('func_value', a.id, ind(k));
                end

            else error('Specify value for x as p(x)');
            end

        end

        function s = get.char(f)
            % Get the formatted string to display the function.
            %
            % s = f.char
            %
            % :param f:
            %     Instance of class :mat:func:`Func`
            % :return:
            %     Formatted string displaying the function
            %
            if strcmp(f.typ, 'sum')
                s = ['(' char(f.f1) ') + (' char(f.f2) ')'];
            elseif strcmp(f.typ, 'diff')
                s = ['(' char(f.f1) ') - (' char(f.f2) ')'];
            elseif strcmp(f.typ, 'prod')
                s = ['(' char(f.f1) ') * (' char(f.f2) ')'];
            elseif strcmp(f.typ, 'ratio')
                s = ['(' char(f.f1) ') / (' char(f.f2) ')'];
            elseif all(p.coeffs == 0)
                s = '0';
            elseif strcmp(f.typ, 'polynomial')
                d = length(p.coeffs) - 1;
                s = [];
                nn = 0;

                for b = f.coeffs;
                    cc(d + 1 - nn) = b;
                    nn = nn + 1;
                end

                for a = cc;

                    if a ~= 0;

                        if ~isempty(s)

                            if a > 0
                                s = [s ' + '];
                            else
                                s = [s ' - '];
                                a = -a;
                            end

                        end

                        if a ~= 1 || d == 0
                            s = [s num2str(a)];

                            if d > 0
                                s = [s '*'];
                            end

                        end

                        if d >= 2
                            s = [s 'x^' int2str(d)];
                        elseif d == 1
                            s = [s 'x'];
                        end

                    end

                    d = d - 1;
                end
            elseif strcmp(f.typ, 'gaussian')
                s = ['Gaussian(' num2str(f.coeffs(1)) ',' ...
                    num2str(f.coeffs(2)) ',' num2str(f.coeffs(3)) ')'];
            elseif strcmp(f.typ, 'fourier')
                c = reshape(f.coeffs, [], 2);
                Ao = c(1, 1);
                w = c(1, 2);
                A = c(2:end, 1);
                B = c(2:end, 2);
                N = size(c, 1) - 1;

                if Ao ~= 0
                    s = num2str(Ao / 2);
                else
                    s = '';
                end

                for n = 1:N

                    if A(n) ~= 0

                        if A(n) < 0
                            prefix = ' - ';
                        elseif s
                            prefix = ' + ';
                        else
                            prefix = '';
                        end

                        s = [s prefix num2str(abs(A(n))), ...
                            '*cos(' num2str(n * w) '*x)'];
                    end

                    if B(n) ~= 0

                        if B(n) < 0
                            prefix = ' - ';
                        elseif s
                            prefix = ' + ';
                        else
                            prefix = '';
                        end

                        s = [s prefix num2str(abs(B(n))), ...
                            '*sin(' num2str(n * w) '*x)'];
                    end

                end
            else
                s = ['*** char not yet implemented for' f.typ ' ***'];
            end

        end

    end

end

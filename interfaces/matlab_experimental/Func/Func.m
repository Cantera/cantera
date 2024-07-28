classdef Func < handle

    properties (SetAccess = immutable)
        id
    end

    methods
        %% Func Class Constructor

        function x = Func(typ, varargin)
            % Func Class ::
            %
            %     >> x = Func(typ, coeff)
            %
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
            %     * Basic functors: ``'sin'``, ``'cos'``, ``'exp'``, ``'log'``,
            %       ``'pow'``, ``'constant'``. Use no or scalar parameters, for example
            %
            %       >> x = Func('cos')
            %       >> x = Func('sin', 2)
            %
            %     * Advanced functors: ``'polynomial'``, ``'Fourier'``, ``'Gaussian'``,
            %       ``'Arrhenius'``. Use vector parameter, for example
            %
            %       >> x = Func('polynomial', [1 2 3])  % x^2 + 2x + 3
            %
            %     * Tabulation functors: ``'tabulated-linear'``,
            %       ``'tabulated-previous'``. Use pair of vector parameters, for example
            %
            %       >> x = Func('tabulated-linear', [0 2 4], [1 3 5])
            %
            %     * Compounding functors: ``'sum'``, ``'diff'``, ``'product'``,
            %       ``'ratio'``, ``'composite'``. Use two functor parameters or a
            %       corresponding operator, for example
            %
            %       >> x = Func('sum', Func('sin', 2.), Func('cos', 3.))
            %       >> x = Func('sin', 2.) + Func('cos', 3.)  % alternative
            %
            %     * Modifying functors: ``'times-constant'``, ``'plus-constant'``,
            %       ``'periodic'``. Use one functor and one scalar, for example
            %
            %       >> x = Func('times-constant', Func('sin', 2.), 2.)
            %
            % Note: this MATLAB class shadows underlying C++ Cantera classes derived
            % from "Func1". See the Cantera C++ documentation for more details.
            %
            % :param typ:
            %     String indicating type of functor to create.
            % :param varargin:
            %     Parameters required for the functor; values depend on definition.
            % :return:
            %     Instance of class :mat:class:`Func`.

            ctIsLoaded;

            if ~isa(typ, 'char')
                error('Function type must be a string');
            end

            if length(varargin) == 0
                % simple functor with no parameter
                x.id = ctFunc('func_new_basic', typ, 1.);
            elseif length(varargin) == 1
                coeffs = varargin{1};
                if length(coeffs) == 1
                    % simple functor with scalar parameter
                    x.id = ctFunc('func_new_basic', typ, coeffs);
                else
                    % advanced functor with array and no parameter
                    x.id = ctFunc('func_new_advanced', typ, length(coeffs), coeffs);
                end
            elseif length(varargin) == 2
                arg1 = varargin{1};
                arg2 = varargin{2};
                if isa(arg1, 'Func') && isa(arg2, 'Func')
                    % compound functor
                    x.id = ctFunc('func_new_compound', typ, arg1.id, arg2.id);
                elseif isa(arg1, 'Func') && isa(arg2, 'double') && length(arg2) == 1
                    % modified functor
                    x.id = ctFunc('func_new_modified', typ, arg1.id, arg2);
                elseif isa(arg1, 'double') && isa(arg2, 'double')
                    % tabulating functors
                    coeffs = [varargin{1}, varargin{2}];
                    x.id = ctFunc('func_new_advanced', typ, length(coeffs), coeffs);
                else
                    error('Invalid arguments');
                end
            else
                error('Invalid arguments');
            end
        end

        %% Func Class Destructor

        function delete(f)
            % Delete the :mat:class:`Func` object.

            ctFunc('func_del', f.id);
        end

        %% Func Class Utility Methods

        function r = plus(f1,f2)
            if isa(f1, 'Func') && isa(f2, 'Func')
                r = Func('sum', f1, f2);
            elseif isa(f1, 'Func') && isa(f2, 'double') && length(f2) == 1
                r = Func('plus-constant', f1, f2)
            else
                error('Invalid arguments')
            end
        end

        function r = minus(f1,f2)
            if isa(f1, 'Func') && isa(f2, 'Func')
                r = Func('diff', f1, f2);
            elseif isa(f1, 'Func') && isa(f2, 'double') && length(f2) == 1
                r = Func('plus-constant', f1, -f2)
            else
                error('Invalid arguments')
            end
        end

        function r = mtimes(f1,f2)
            if isa(f1, 'Func') && isa(f2, 'Func')
                r = Func('product', f1, f2);
            elseif isa(f1, 'Func') && isa(f2, 'double') && length(f2) == 1
                r = Func('times-constant', f1, f2)
            else
                error('Invalid arguments')
            end
        end

        function r = mrdivide(f1,f2)
            r = Func('ratio', f1, f2);
        end

        function s = type(f)
            % Return function type.
            s = ctString('func_type', f.id);
        end

        function display(f)
            % Display the equation of the input function on the terminal.

            disp(' ');
            disp([inputname(1), ' = '])
            disp(' ');
            disp(['   ' f.write])
            disp(' ');
        end

        function b = subsref(a, s)
            % Redefine subscripted references for functors. ::
            %
            %     >> b = a.subsref(s)
            %
            % :param a:
            %     Instance of class :mat:class:`Func`.
            % :param s:
            %     Value at which the function should be evaluated.
            % :return:
            %     Returns the value of the function evaluated at ``s``.

            if length(s) > 1
                aa = eval(['a.', s(1).subs]);
                b = subsref(aa, s(2:end));
                return
            end
            if strcmp(s.type, '()')
                ind = s.subs{:};
                b = zeros(1, length(ind));

                for k = 1:length(ind)
                    b(k) = ctFunc('func_value', a.id, ind(k));
                end

            elseif strcmp(s.type, '.')
                b = eval(['a.', s.subs]);
            else error('Specify value for x as p(x)');
            end

        end

        function s = write(f)
            % Get the formatted string to display the function. ::
            %
            %     >> s = f.write
            %
            % :param f:
            %     Instance of class :mat:class:`Func`.
            % :return:
            %     LaTeX-formatted string displaying the function.
            arg = 'x';
            s = ctString('func_write3', f.id, arg);
        end

    end

end

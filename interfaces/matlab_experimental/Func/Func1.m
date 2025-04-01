classdef Func1 < handle

    properties (SetAccess = immutable)
        id
    end

    properties (SetAccess = protected)
        type
    end

    methods
        %% Func1 Class Constructor

        function x = Func1(typ, varargin)
            % Func1 Class ::
            %
            %     >> x = Func1(typ, coeff)
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
            %       >> x = Func1('cos')
            %       >> x = Func1('sin', 2)
            %
            %     * Advanced functors: ``'polynomial3'``, ``'Fourier'``, ``'Gaussian'``,
            %       ``'Arrhenius'``. Use vector parameter, for example
            %
            %       >> x = Func1('polynomial3', [1 2 3])  % x^2 + 2x + 3
            %
            %     * Tabulation functors: ``'tabulated-linear'``,
            %       ``'tabulated-previous'``. Use pair of vector parameters, for example
            %
            %       >> x = Func1('tabulated-linear', [0 2 4], [1 3 5])
            %
            %     * Compounding functors: ``'sum'``, ``'diff'``, ``'product'``,
            %       ``'ratio'``, ``'composite'``. Use two functor parameters or a
            %       corresponding operator, for example
            %
            %       >> x = Func1('sum', Func1('sin', 2.), Func1('cos', 3.))
            %       >> x = Func1('sin', 2.) + Func1('cos', 3.)  % alternative
            %
            %     * Modifying functors: ``'times-constant'``, ``'plus-constant'``,
            %       ``'periodic'``. Use one functor and one scalar, for example
            %
            %       >> x = Func1('times-constant', Func1('sin', 2.), 2.)
            %
            % Note: this MATLAB class shadows underlying C++ Cantera classes derived
            % from "Func1". See the Cantera C++ documentation for more details.
            %
            % :param typ:
            %     String indicating type of functor to create.
            % :param varargin:
            %     Parameters required for the functor; values depend on definition.
            % :return:
            %     Instance of class :mat:class:`Func1`.

            ctIsLoaded;

            if isnumeric(typ)
                % instantiate from handle
                x.id = typ;
                return
            end

            if ~isa(typ, 'char')
                error('Function type must be a string');
            end
            func1Type = ctString('func_check', typ);

            if func1Type == "undefined"
                error(['Functor ''' typ ''' is not implemented'])
            end

            if length(varargin) == 0 && func1Type == "standard"
                % basic functor with no parameter
                x.id = ctFunc('func_new_basic', typ, 1.);
            elseif length(varargin) == 1 && func1Type == "standard"
                coeffs = varargin{1};
                if length(coeffs) == 1
                    % basic functor with scalar parameter
                    x.id = ctFunc('func_new_basic', typ, coeffs);
                elseif isa(coeffs, 'double')
                    % advanced functor with array parameter
                    x.id = ctFunc('func_new_advanced', typ, length(coeffs), coeffs);
                else
                    error('Invalid arguments for functor')
                end
            elseif length(varargin) == 2
                arg1 = varargin{1};
                arg2 = varargin{2};
                if func1Type == "compound"
                    % compounding functor
                    if isa(arg1, 'double') && length(arg1) == 1
                        arg1 = Func1('constant', arg1);
                    elseif isa(arg2, 'double') && length(arg2) == 1
                        arg2 = Func1('constant', arg2);
                    end
                    if ~isa(arg1, 'Func1') || ~isa(arg2, 'Func1')
                        error('Invalid arguments for compounding functor')
                    end
                    x.id = ctFunc('func_new_compound', typ, arg1.id, arg2.id);
                elseif func1Type == "modified"
                    % modifying functor
                    if ~isa(arg1, 'Func1') || ~isa(arg2, 'double') || length(arg2) > 1
                        error('Invalid arguments for modifying functor')
                    end
                    x.id = ctFunc('func_new_modified', typ, arg1.id, arg2);
                else % func1Type == "standard"
                    % tabulating functors
                    if ~isa(arg1, 'double') || ~isa(arg2, 'double')
                        error('Invalid arguments for tabulating functor')
                    end
                    coeffs = [varargin{1}, varargin{2}];
                    x.id = ctFunc('func_new_advanced', typ, length(coeffs), coeffs);
                end
            else
                error('Invalid number of arguments');
            end
        end

        %% Func1 Class Destructor

        function delete(f)
            % Delete the :mat:class:`Func1` object.

            if isempty(f.id)
                return
            end
            ctFunc('func_del', f.id);
        end

        %% Func1 Class Utility Methods

        function r = plus(f1,f2)
            if isa(f1, 'double') && length(f1) == 1
                f1 = Func1('constant', f1);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = Func1('constant', f2);
            end
            id = ctFunc('func_new_sum', f1.id, f2.id);
            r = Func1(id);
        end

        function r = minus(f1,f2)
            if isa(f1, 'double') && length(f1) == 1
                f1 = Func1('constant', f1);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = Func1('constant', f2);
            end
            id = ctFunc('func_new_diff', f1.id, f2.id);
            r = Func1(id);
        end

        function r = mtimes(f1,f2)
            if isa(f1, 'double') && length(f1) == 1
                f1 = Func1('constant', f1);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = Func1('constant', f2);
            end
            id = ctFunc('func_new_prod', f1.id, f2.id);
            r = Func1(id);
        end

        function r = mrdivide(f1,f2)
            if isa(f1, 'double') && length(f1) == 1
                f1 = Func1('constant', f1);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = Func1('constant', f2);
            end
            id = ctFunc('func_new_ratio', f1.id, f2.id);
            r = Func1(id);
        end

        function s = get.type(f)
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
            %     Instance of class :mat:class:`Func1`.
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
            %     Instance of class :mat:class:`Func1`.
            % :return:
            %     LaTeX-formatted string displaying the function.
            arg = 'x';
            s = ctString('func_write', f.id, arg);
        end

    end

end

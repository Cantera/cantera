classdef Func1 < handle

    properties (SetAccess = immutable)
        id = -1
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

            ct.isLoaded(true);

            if isnumeric(typ)
                % instantiate from handle
                x.id = typ;
                return
            end

            if ~isa(typ, 'char')
                error('Function type must be a string');
            end
            func1Type = ct.impl.getString('mFunc1_checkFunc1', typ);

            if func1Type == "undefined"
                error(['Functor ''' typ ''' is not implemented'])
            end

            if length(varargin) == 0 && func1Type == "standard"
                % basic functor with no parameter
                x.id = ct.impl.call('mFunc1_newBasic', typ, 1.);
            elseif length(varargin) == 1 && func1Type == "standard"
                coeffs = varargin{1};
                if length(coeffs) == 1
                    % basic functor with scalar parameter
                    x.id = ct.impl.call('mFunc1_newBasic', typ, coeffs);
                elseif isa(coeffs, 'double')
                    % advanced functor with array parameter
                    x.id = ct.impl.call('mFunc1_newAdvanced', typ, coeffs);
                else
                    error('Invalid arguments for functor')
                end
            elseif length(varargin) == 2
                arg1 = varargin{1};
                arg2 = varargin{2};
                if func1Type == "compound"
                    % compounding functor
                    if isa(arg1, 'double') && length(arg1) == 1
                        arg1 = ct.Func1('constant', arg1);
                    elseif isa(arg2, 'double') && length(arg2) == 1
                        arg2 = ct.Func1('constant', arg2);
                    end
                    if ~isa(arg1, 'ct.Func1') || ~isa(arg2, 'ct.Func1')
                        error('Invalid arguments for compounding functor')
                    end
                    x.id = ct.impl.call('mFunc1_newCompound', typ, arg1.id, arg2.id);
                elseif func1Type == "modified"
                    % modifying functor
                    if ~isa(arg1, 'ct.Func1') || ~isa(arg2, 'double') || length(arg2) > 1
                        error('Invalid arguments for modifying functor')
                    end
                    x.id = ct.impl.call('mFunc1_newModified', typ, arg1.id, arg2);
                else % func1Type == "standard"
                    % tabulating functors
                    if ~isa(arg1, 'double') || ~isa(arg2, 'double')
                        error('Invalid arguments for tabulating functor')
                    end
                    coeffs = [varargin{1}, varargin{2}];
                    x.id = ct.impl.call('mFunc1_newAdvanced', typ, coeffs);
                end
            else
                error('Invalid number of arguments');
            end
        end

        %% Func1 Class Destructor

        function delete(obj)
            % Delete the :mat:class:`Func1` object.
            if obj.id >= 0
                ct.impl.call('mFunc1_del', obj.id);
            end
        end

        %% Func1 Class Utility Methods

        function r = plus(obj,f2)
            if isa(obj, 'double') && length(obj) == 1
                obj = ct.Func1('constant', obj);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = ct.Func1('constant', f2);
            end
            id = ct.impl.call('mFunc1_newSumFunction', obj.id, f2.id);
            r = ct.Func1(id);
        end

        function r = minus(obj,f2)
            if isa(obj, 'double') && length(obj) == 1
                obj = ct.Func1('constant', obj);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = ct.Func1('constant', f2);
            end
            id = ct.impl.call('mFunc1_newDiffFunction', obj.id, f2.id);
            r = ct.Func1(id);
        end

        function r = mtimes(obj,f2)
            if isa(obj, 'double') && length(obj) == 1
                obj = ct.Func1('constant', obj);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = ct.Func1('constant', f2);
            end
            id = ct.impl.call('mFunc1_newProdFunction', obj.id, f2.id);
            r = ct.Func1(id);
        end

        function r = mrdivide(obj,f2)
            if isa(obj, 'double') && length(obj) == 1
                obj = ct.Func1('constant', obj);
            elseif isa(f2, 'double') && length(f2) == 1
                f2 = ct.Func1('constant', f2);
            end
            id = ct.impl.call('mFunc1_newRatioFunction', obj.id, f2.id);
            r = ct.Func1(id);
        end

        function s = get.type(obj)
            % Return function type.
            s = ct.impl.getString('mFunc1_type', obj.id);
        end

        function display(obj)
            % Display the equation of the input function on the terminal.

            disp(' ');
            disp([inputname(1), ' = '])
            disp(' ');
            disp(['   ' obj.write])
            disp(' ');
        end

        function b = subsref(obj, s)
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
                name = s(1).subs;
                aa = obj.(name);
                b = subsref(aa, s(2:end));
                return
            end
            if strcmp(s.type, '()')
                ind = s.subs{:};
                b = zeros(1, length(ind));

                for k = 1:length(ind)
                    b(k) = ct.impl.call('mFunc1_eval', obj.id, ind(k));
                end

            elseif strcmp(s.type, '.')
                name = s.subs;
                b = obj.(name);
            else
                error('Specify value for x as p(x)');
            end

        end

        function s = write(obj)
            % Get the formatted string to display the function. ::
            %
            %     >> s = f.write
            %
            % :param f:
            %     Instance of class :mat:class:`Func1`.
            % :return:
            %     LaTeX-formatted string displaying the function.
            arg = 'x';
            s = ct.impl.getString('mFunc1_write', obj.id, arg);
        end

    end

end

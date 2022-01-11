classdef Func < handle

    properties
        f1
        f2
        coeffs
        id
        typ
    end

    methods
        %% Functor class constructor

        function x = Func(typ, n, p)
            % Functor class constructor.
            %
            % A functor is an object that behaves like a function. Cantera
            % defines a set of functors to use to create arbitrary
            % functions to specify things like heat fluxes, piston speeds,
            % etc., in reactor network simulations.
            %
            % The main feature of a functor class is that it overloads the
            % '()' operator to evaluate the function. For example, suppose
            % object 'f' is a functor that evaluates the polynomial
            % :math:'2x^2 - 3x + 1'. Then writing 'f(2)' would cause the
            % method that evaluates the function to be invoked, and would
            % pass it the argument '2'. The return value would be 3.
            %
            % The types of functors you can create in Cantera are these:
            % 1. A polynomial;
            % 2. A Fourier series;
            % 3. A sum of Arrhenius terms;
            % 4. A Gaussian.
            %
            % You can also create composite functors by adding,
            % multiplying, or dividing these basic functors or other
            % composite functors.
            %
            % Note: this MATLAB class shadows the underlying C++ Cantera
            % class 'Func1'. See the Cantera C++ documentation for more
            % details.
            %
            % parameter typ:
            %    String indicating type of functor to create. Possible
            %    values are:
            %    * 'polynomial'
            %    * 'fourier'
            %    * 'gaussian'
            %    * 'arrhenius'
            %    * 'sum'
            %    * 'diff'
            %    * 'ratio'
            %    * 'composite'
            %    * 'periodic'
            %
            % parameter n:
            %    Number of parameters required for the functor
            % parameter p:
            %    Vector of parameters
            % return:
            %    Instance of class :mat:func:`Func`

            checklib;

            if ~isa(typ, 'char')
                error('Function type must be a string');
            end

            x.f1 = 0;
            x.f2 = 0;
            x.coeffs = 0;

            function nn = new_func(itype, n, p)
                % helper function to pass the correct parameters to the C
                % library
                if itype < 20
                    ptr = libpointer('doublePtr', p);
                    [m, n] = size(p);
                    lenp = m * n;
                    nn = calllib(ct, 'func_new', type, n, lenp, ptr);
                elseif itype < 45
                    m = n;
                    nn = calllib(ct, 'func_new', type, n, m, 0);
                else
                    ptr = libpointer('doublePtr', p);
                    nn = calllib(ct, 'func_new', type, n, 0, ptr);
                end
            end

            if itype > 0
                x.coeffs = p;
                x.id = new_func(itype, n, p);
            elseif strcmp(typ, 'periodic')
                itype = 50;
                x.f1 = n;
                x.coeffs = p;
                x.id = new_func(itype, n.id, p);
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
                x.id = new_func(itype, n.id, p.id);
            end

            x.typ = typ;
        end

        %% Utility methods

        function clear(f)
            % Clear the functor from memory.
            checklib;
            calllib(ct, 'func_del', f.id);
        end

        function display(a)
            % Display the equation of the input function on the terminal.
            %
            % parameter a:
            %    Instance of class 'Func'

            disp(' ');
            disp([inputname(1),' = '])
            disp(' ');
            disp(['   ' char(a)])
            disp(' ');
        end

        %% Functor methods

        function r = plus(a, b)
            % Get a functor representing the sum two input functors 'a' and
            % 'b'.
            r = Func('sum', a, b);
        end

        function r = rdivide(a, b)
            % Get a functor that is the ratio of the input functors 'a' and
            % 'b'.
            r = Func('ratio', a, b);
        end

        function r = times(a, b)
             % Get a functor that multiplies two functors 'a' and 'b'
             r = Func('prod', a, b);
        end

        function b = subsref(a, s)
            % Redefine subscripted references for functors.
            %
            % parameter a:
            %    Instance of class 'Func'
            % parameter s:
            %    Value at which the function should be evaluated.
            % return:
            %    The value of the function evaluated at 's'.

            checklib;
            if strcmp(s.type, '()')
                ind= s.subs{:};
                 b = zeros(1, length(ind));
                 for k = 1:length(ind)
                     b(k) = calllib(ct, 'func_value', a.id, ind(k));
                 end
            else error('Specify value for x as p(x)');
            end
        end

        function s = char(p)
           % Get the formatted string to display the function.
           if strcmp(p.typ,'sum')
               s = ['(' char(p.f1) ') + (' char(p.f2) ')'];
           elseif strcmp(p.typ,'diff')
               s = ['(' char(p.f1) ') - (' char(p.f2) ')'];
           elseif strcmp(p.typ,'prod')
               s = ['(' char(p.f1) ') * (' char(p.f2) ')'];
           elseif strcmp(p.typ,'ratio')
               s = ['(' char(p.f1) ') / (' char(p.f2) ')'];
           elseif all(p.coeffs == 0)
               s = '0';
           else
               if strcmp(p.typ,'polynomial')
                   d = length(p.coeffs) - 1;
                   s = [];
                   nn = 0;
                   for b = p.coeffs;
                       cc(d+1-nn) = b;
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
               elseif strcmp(p.typ, 'gaussian')
                   s = ['Gaussian(' num2str(p.coeffs(1)) ',' ...
                       num2str(p.coeffs(2)) ',' ...
                       num2str(p.coeffs(3)) ')'];
               elseif strcmp(p.typ, 'fourier')
                   c = reshape(p.coeffs, [], 2);
                   Ao = c(1, 1);
                   w = c(1, 2);
                   A = c(2:end, 1);
                   B = c(2:end, 2);
                   N = size(c, 1) - 1;
                   if Ao ~= 0
                       s = num2str(Ao/2);
                   else
                       s = '';
                   end
                   for n=1:N
                       if A(n) ~= 0
                           if A(n) < 0
                               prefix = ' - ';
                           elseif s
                               prefix = ' + ';
                           else
                               prefix = '';
                           end
                       s = [s prefix num2str(abs(A(n))), ...
                            '*cos(' num2str(n*w) '*x)'];
                       end
                       if B(n) ~= 0
                           if B(n) < 0
                               prefix = ' - ';
                           elseif s
                               prefix = ' + ';
                           else
                               prefix = '';
                           end
                       s = [s prefix num2str(abs(B(n))),...
                            '*sin(' num2str(n*w) '*x)'];
                       end
                   end
               else
                   s = ['*** char not yet implemented for' p.typ ' ***'];
               end
           end
        end

    end
end

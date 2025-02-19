classdef ctTestFunc1 < matlab.unittest.TestCase

    properties (SetAccess = immutable)
        rtol = 1e-6;
        atol = 1e-8;
    end

    methods (TestClassSetup)

        function testSetUp(self)
            ctTestSetUp
        end

    end

    methods (TestClassTeardown)

        function testTearDown(self)
            ctCleanUp
            ctTestTearDown
        end

    end

    methods (Test)

        function testFunction(self)
            f = Func1('sin', 1);
            for i = [0, 0.1, 0.7]
                self.verifyEqual(f(i), sin(i), 'AbsTol', self.atol);
            end

            clear f
        end

        function testLambda(self)
            f1 = Func1('sin', 1);
            f2 = Func1('pow', 0.5);
            f = Func1('product', f1, f2);

            for i = [0.1, 0.7, 4.5]
                self.verifyEqual(f(i), sin(i)*sqrt(i), 'AbsTol', self.atol);
            end

            clear f1 f2 f
        end

        function testConstant(self)
            f = Func1('constant', 5);

            for i = [0.1, 0.7, 4.5]
                self.verifyEqual(f(i), 5, 'AbsTol', self.atol);
            end

            clear f
        end

        function testSimple(self)
            functors = {'sin', 'cos', 'exp', 'log'};
            coeff = 2.34;

            for name = functors
                f = Func1(name{:}, coeff);
                self.verifySubstring(f.write, name{:});
                self.verifyEqual(f.type, name{:});
                for val = [0.1, 1, 10]
                    self.verifyEqual(f(val), feval(name{:}, coeff*val), ...
                                     'AbsTol', self.atol);
                end
                clear f
            end
        end

        function testCompound(self)
            f0 = 3.1415;
            f1 = Func1('pow', 2);
            f2 = Func1('sin');
            val = [0.1, 1, 10];
            functors = containers.Map(...
                {'sum', 'diff', 'product', 'ratio'}, ...
                {@(x, y) x + y, @(x, y) x - y, @(x, y) x * y, @(x, y) x / y});

            for k = keys(functors)
                func = Func1(k{:}, f1, f2);
                f = functors(k{:});
                self.verifyFalse(contains(func.write, k{:}));
                self.verifyEqual(k{:}, func.type);
                for v = val
                    x1 = func(v);
                    x2 = f(f1(v), f2(v));
                    self.verifyEqual(x1, x2, 'absTol', self.atol);
                end
                clear func
            end

            f3 = Func1('sum', f0, f1);
            f4 = Func1('sum', f2, f0);
            for v = val
                self.verifyEqual(f3(v), f0 + f1(v), 'absTol', self.atol);
                self.verifyEqual(f4(v), f0 + f2(v), 'absTol', self.atol);
            end

            try
                f5 = Func1('sum', f0, f0);
            catch ME
                self.verifySubstring(ME.message, 'Invalid arguments');
            end

            try
                f6 = Func1('sum', 'spam', 'eggs');
            catch ME
                self.verifySubstring(ME.message, 'Invalid arguments');
            end

            clear f1 f2 f3 f4
        end

        function testNewSum(self)
            f1 = Func1('sin', 3);

            f2 = f1 + 0;
            self.verifyEqual(f2.type, 'sin');
            self.verifyEqual(f2(0.5), sin(0.5 * 3), 'absTol', self.atol);
            clear f2

            f2 = 0 + f1;
            self.verifyEqual(f2.type, 'sin');
            self.verifyEqual(f2(0.5), sin(0.5 * 3), 'absTol', self.atol);
            clear f2

            f2 = f1 + Func1('constant', 1);
            self.verifyEqual(f2.type, 'plus-constant');
            self.verifyEqual(f2(0.5), sin(0.5 * 3) + 1, 'absTol', self.atol);
            clear f2

            f2 = f1 + f1;
            self.verifyEqual(f2.type, 'times-constant');
            self.verifyEqual(f2(0.5), sin(0.5 * 3) * 2, 'absTol', self.atol);
            clear f2

            clear f1
        end

        function testNewDiff(self)
            f1 = Func1('sin', 3);

            f2 = f1 - 0;
            self.verifyEqual(f2.type, 'sin');
            self.verifyEqual(f2(0.5), sin(0.5 * 3), 'absTol', self.atol);
            clear f2

            f2 = 0 - f1;
            self.verifyEqual(f2.type, 'times-constant');
            self.verifyEqual(f2(0.5), -sin(0.5 * 3), 'absTol', self.atol);
            clear f2

            f2 = f1 - Func1('constant', 1);
            self.verifyEqual(f2.type, 'plus-constant');
            self.verifyEqual(f2(0.5), sin(0.5 * 3) - 1, 'absTol', self.atol);
            clear f2

            f2 = f1 - f1;
            self.verifyEqual(f2.type, 'constant');
            self.verifyEqual(f2(0.5), 0, 'absTol', self.atol);
            clear f2

            clear f1
        end

        function testNewProd(self)
            f1 = Func1('sin', 3);

            f2 = f1 * Func1('constant', 1);
            self.verifyEqual(f2.type, 'sin');
            self.verifyEqual(f2(0.5), sin(0.5 * 3), 'absTol', self.atol);
            clear f2

            f2 = f1 * 2;
            self.verifyEqual(f2.type, 'times-constant');
            self.verifyEqual(f2(0.5), sin(0.5 * 3) * 2, 'absTol', self.atol);
            clear f2

            f2 = f1 * 0;
            self.verifyEqual(f2.type, 'constant');
            self.verifyEqual(f2(0.5), 0, 'absTol', self.atol);
            clear f2

            clear f1
        end

        function testNewRatio(self)
            f1 = Func1('sin', 3);

            f2 = f1 / Func1('constant', 1);
            self.verifyEqual(f2.type, 'sin');
            self.verifyEqual(f2(0.5), sin(0.5 * 3), 'absTol', self.atol);
            clear f2

            f2 = 0 / f1;
            self.verifyEqual(f2.type, 'constant');
            self.verifyEqual(f2(0.5), 0, 'absTol', self.atol);
            clear f2

            f2 = f1 / 2;
            self.verifyEqual(f2.type, 'times-constant');
            self.verifyEqual(f2(0.5), sin(0.5 * 3) / 2, 'absTol', self.atol);
            clear f2

            try
                f2 = f1 / 0;
            catch ME
                self.verifySubstring(ME.identifier, 'Cantera:ctError');
                self.verifySubstring(ME.message, 'Division by zero');
                clear f2
            end

            clear f1
        end

        function testModified(self)
            f1 = Func1('sin');
            const = 2.34;
            val = [0.1, 1, 10];
            functors = containers.Map(...
                {'plus-constant', 'times-constant'}, ...
                {@(x, y) x + y, @(x, y) x * y});

            for k = keys(functors)
                func = Func1(k{:}, f1, const);
                f = functors(k{:});
                self.verifyFalse(contains(func.write, k{:}));
                self.verifyEqual(k{:}, func.type);
                for v = val
                    x1 = func(v);
                    x2 = f(f1(v), const);
                    self.verifyEqual(x1, x2, 'absTol', self.atol);
                end
                clear func

                try
                    func = Func1(k{:}, const, f1);
                catch ME
                    self.verifySubstring(ME.message, 'Invalid arguments');
                end
            end
        end

        function testTabulated(self)
            t = [0, 1, 2];
            fval = [2, 1, 0];
            v1 = [-0.5, 0, 0.5, 1.0, 1.5, 2, 2.5];
            v2 = [2.0, 2.0, 1.5, 1.0, 0.5, 0.0, 0.0];
            v3 = [2.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0];

            func1 = Func1('tabulated-linear', t, fval);
            self.verifyEqual(func1.type, 'tabulated-linear');
            for i = length(v1)
                self.verifyEqual(func1(v1(i)), v2(i), 'absTol', self.atol);
            end

            func2 = Func1('tabulated-previous', t, fval);
            self.verifyEqual(func2.type, 'tabulated-previous');
            for i = length(v1)
                self.verifyEqual(func2(v1(i)), v3(i), 'absTol', self.atol);
            end

            try
                func3 = Func1('tabulated-linear', 0:2, 0:3);
            catch ME
                self.verifySubstring(ME.message, 'even number of entries');
            end

            try
                func3 = Func1('tabulated-linear', [], []);
            catch ME
                self.verifySubstring(ME.message, 'at least 4 entries');
            end

            try
                func3 = Func1('tabulated-linear', [0, 1, 0.5, 2], [2, 1, 1, 0]);
            catch ME
                self.verifySubstring(ME.message, 'monotonically');
            end

            clear func1 func2
        end

    end

end

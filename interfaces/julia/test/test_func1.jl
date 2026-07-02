using Cantera
using Test

# Reference values cross-checked against Python `cantera.Func1` (v3.2.0).

const CT = Cantera

@testset "Func1" begin
    @testset "basic evaluation" begin
        f = CT.Func1("sin", 2.0)
        @test f(pi / 4) ≈ sin(2 * pi / 4)
        @test CT.evaluate(f, pi / 4) ≈ sin(2 * pi / 4)
        @test CT.func_type(f) == "sin"

        c = CT.constant_function(3.5)
        @test c(123.0) ≈ 3.5
    end

    @testset "advanced (polynomial)" begin
        # polynomial3 with [1,2,3,4] -> t^3 + 2t^2 + 3t + 4
        p = CT.Func1("polynomial3", [1.0, 2.0, 3.0, 4.0])
        @test p(2.0) ≈ 26.0     # Python: 26.0
        @test p(0.0) ≈ 4.0
    end

    @testset "arithmetic" begin
        s = CT.Func1("constant", 2.0) + CT.Func1("constant", 3.0)
        @test s(1.0) ≈ 5.0
        d = CT.Func1("constant", 5.0) - CT.Func1("constant", 2.0)
        @test d(1.0) ≈ 3.0
        pr = CT.Func1("constant", 2.0) * CT.Func1("constant", 4.0)
        @test pr(1.0) ≈ 8.0
        r = CT.Func1("constant", 6.0) / CT.Func1("constant", 3.0)
        @test r(1.0) ≈ 2.0
    end

    @testset "derivative" begin
        # d/dt sin(w t) = w cos(w t); at t=0 -> w
        f = CT.Func1("sin", 2.0)
        df = CT.derivative(f)
        @test df(0.0) ≈ 2.0
    end

    @testset "show / write" begin
        f = CT.Func1("sin", 2.0)
        @test CT.write(f) isa String
        @test occursin("sin", sprint(show, f))
    end
end

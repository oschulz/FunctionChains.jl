# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test

using AffineMaps
using InverseFunctions
using ChangesOfVariables

include("testfuncs.jl")

@testset "as_function" begin
    mc = MulCallable(2.7)
    mf = Mul(2.9)
    tc = Complex

    x = 4.2

    @test @inferred(asfunction(log)) === log
    @test @inferred(asfunction(mc)) isa Function

    @test @inferred(fbcast(asfunction(Int))) === fbcast(Int)
    @test @inferred(inverse(asfunction(Int))) === inverse(Int)
    @test @inferred(with_logabsdet_jacobian(asfunction(Int), 7)) === with_logabsdet_jacobian(Int, 7)


    for f in [mc, mf, tc]
        @test @inferred(asfunction(f)) isa Function
        if f isa Function
            @test @inferred(asfunction(f)) === f
        else
            @test @inferred(asfunction(f)) isa Function
        end

        if !(f isa Type)
            @test @inferred(typed_callable(f)) === f
        else
            @test @inferred(typed_callable(f)) isa Function
        end

        ff = asfunction(f)
        @test @inferred(ff(x)) == f(x)

        tc = typed_callable(f)
        @test @inferred(tc(x)) == f(x)
    end

    @test @inferred(asfunction(mc)) isa FunctionChains.AsFunction
    f_mc = asfunction(mc)
    y = 5.2
    @test @inferred(f_mc(x)) == mc(x)
    @test @inferred(f_mc(x, y)) == mc(x, y)
    @test @inferred(f_mc(x, k = 2)) == mc(x, k = 2)
    @test @inferred(f_mc(x, y, k = 2)) == mc(x, y, k = 2)

    asf_mul = FunctionChains.AsFunction(Mul(rand(2,2)))
    @test asf_mul isa FunctionChains.AsFunction
    asf_mul_2 = deepcopy(asf_mul)
    @test !(asf_mul.f.A === asf_mul_2.f.A)
    @test @inferred(asf_mul == asf_mul_2)
    @test @inferred(isapprox(asf_mul, asf_mul_2, rtol = 1e-5))
end

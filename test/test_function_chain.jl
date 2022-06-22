# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test
using Base.Iterators: repeated, take
using InverseFunctions, ChangesOfVariables

include("getjacobian.jl")
include("testfuncs.jl")


@testset "function_chain" begin
    function test_function_chain(fc, x, invertible::Bool, with_ladj::Bool, label::AbstractString)
        @testset "$label" begin
            @test !any(x -> x isa ComposedFunction, fc.fs)

            cf = foldl(∘, reverse(collect(fc.fs)))
            @test @inferred(fc(x)) == cf(x)

            y_out = @inferred fc(x)
            ys_out = @inferred with_intermediate_results(fc, x)
        
            y_ref = x   
            for (f_i, y_out_i) in zip(fc.fs, ys_out)
                y_ref = f_i(y_ref)
                @test y_out_i == y_ref
            end
            @test y_out == y_ref

            if invertible
                @test @inferred(inverse(fc)) isa FunctionChain
                inv_cf = inverse(cf)
                @test @inferred(inverse(fc)(x)) == inv_cf(x)
                InverseFunctions.test_inverse(fc, x)
            else
                @test @inferred(inverse(fc)) isa NoInverse
            end

            if with_ladj
                y_ladj_ref = with_logabsdet_jacobian(cf, x)
                @test @inferred(with_logabsdet_jacobian(fc, x)) == y_ladj_ref
                ChangesOfVariables.test_with_logabsdet_jacobian(fc, x, getjacobian)
            else
                @test @inferred(with_logabsdet_jacobian(fc, x)) isa NoLogAbsDetJacobian
            end
        end
    end

    @test @inferred(fchain(AffineStep.((*, +, \, -), (3, 2, 2.5, 4)))) isa FunctionChain{<:NTuple{N,AffineStep} where N}
    f = fchain(AffineStep.((*, +, \, -), (3, 2, 2.5, 4)))
    @test @inferred(f([1.1, 2.2])) == 2.5 \ (3 * [1.1, 2.2] .+ 2) .- 4
    test_function_chain(f, [1.1, 2.2], true, true, "Tuple of AffineStep")

    @test @inferred(fchain(AffineStep.(*, 2:5))) isa FunctionChain{Vector{AffineStep{typeof(*),Int}}}
    f = fchain(AffineStep.(*, 2:5))
    @test @inferred(f(1.2)) ≈ 1.2 * prod(2:5)
    test_function_chain(f, 1.2, true, true, "Vector of AffineStep, scalar arg")
    @test @inferred(with_intermediate_results(f, [1.1, 2.2])) ≈ [[2.2, 4.4], [6.6, 13.2], [26.4, 52.8], [132.0, 264.0]]
    test_function_chain(f, [1.1, 2.2], true, true, "Tuple of AffineStep, vector arg")

    @test @inferred(fchain(AffineStep(*, i) for i in 2:5)) isa FunctionChain{<:Base.Generator}
    f = fchain(AffineStep(*, i) for i in 2:5)
    @test @inferred(with_intermediate_results(f, [1.1, 2.2])) ≈ [[2.2, 4.4], [6.6, 13.2], [26.4, 52.8], [132.0, 264.0]]
    test_function_chain(f, [1.1, 2.2], false, true, "Generator of AffineStep")

    @test @inferred(fchain(repeated(AffineStep(*, 2), 5))) isa FunctionChain{<:Base.Iterators.Take}
    f = fchain(repeated(AffineStep(*, 2), 5))
    test_function_chain(f, [1.1, 2.2], true, true, "Generator of AffineStep")

    @test @inferred(inv ∘ fchain() ∘ (expm1 ∘ log)) == fchain((log, expm1, inv))
    @test @inferred(fchain((log, expm1, inv))(0.3)) == (inv ∘ expm1 ∘ log)(0.3)
    test_function_chain(fchain((log, expm1, inv)), 0.3, true, true, "log-expm1-inv")
    test_function_chain(fchain(fill(expm1, 3)), 0.3, true, true, "fill-expm1")
    test_function_chain(fchain((log, sin, expm1, inv)), 0.3, false, false, "log-sin-expm1-inv")
    test_function_chain(fchain(fill(sin, 3)), 0.3, false, false, "fill-sin")
end

# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test
using Base.Iterators: repeated, take
using InverseFunctions, ChangesOfVariables
import Adapt, Functors
using AffineMaps
import ForwardDiff

import FlexiMaps

include("getjacobian.jl")
include("testfuncs.jl")

@testset "function_chain" begin
    fc2cf(fc::FunctionChain) = foldl(∘, reverse(collect(fchainfs(fc))))
    fbcast_fc2cf(fc::FunctionChain) = foldl(∘, reverse(map(Broadcast.BroadcastFunction, collect(fchainfs(fc)))))

    # with_intermediate_results reference implementation:
    ref_wir(fc::FunctionChain{<:Tuple}, x) = foldl((acc, f) -> (acc..., f(last(acc))), values(fchainfs(fc)), init = (x,))[2:end]
    ref_wir(fc::FunctionChain{<:NamedTuple{names}}, x) where names = NamedTuple{names}(ref_wir(fchain(values(fchainfs(fc))...), x))
    ref_wir(fc::FunctionChain, x) = foldl((acc, f) -> [acc..., f(last(acc))], values(fchainfs(fc)), init = [x])[2:end]
    ref_wir_bc(fc::FunctionChain{<:Tuple}, x) = foldl((acc, f) -> (acc..., f.(last(acc))), values(fchainfs(fc)), init = (x,))[2:end]
    ref_wir_bc(fc::FunctionChain{<:NamedTuple{names}}, x) where names = NamedTuple{names}(ref_wir_bc(fchain(values(fchainfs(fc))...), x))
    ref_wir_bc(fc::FunctionChain, x) = foldl((acc, f) -> [acc..., f.(last(acc))], values(fchainfs(fc)), init = [x])[2:end]

    function test_function_chain(fc, x, xs, invertible::Bool, with_ladj::Bool, label::AbstractString)
        @testset "$label" begin
            @test @inferred(fchainfs(fc)) === getfield(fc, :_fs)

            @test !any(x -> x isa ComposedFunction, fchainfs(fc))

            @test length(fc) == length(fchainfs(fc))
            @test (fc...,) == (fchainfs(fc)...,)
            @test [fc...] == [fchainfs(fc)...]

            cf = fc2cf(fc)
            @test @inferred(fc(x)) == cf(x)

            y_out = @inferred fc(x)
            ys_out = @inferred with_intermediate_results(fc, x)

            @test ys_out == ref_wir(fc, x)

            ys_out_values = if fchainfs(fc) isa NamedTuple
                @test propertynames(ys_out) == propertynames(fchainfs(fc))
                values(ys_out)
            else
                ys_out
            end
        
            y_ref = x   
            for (f_i, y_out_i) in zip(fchainfs(fc), ys_out_values)
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

            if !isnothing(xs)
                refbc_f_tpl = fbcast_fc2cf(fc)
                @test @inferred(broadcast(fc, (xs))) == cf.(xs)
                @test @inferred(fbcast(fc)(xs)) == cf.(xs)
                @test @inferred(with_intermediate_results(fbcast(fc), xs)) == ref_wir_bc(fc, xs)
                if with_ladj
                    ly, ladj = @inferred(with_logabsdet_jacobian(fbcast(fc), xs))
                    ly_ref, ladj_ref = with_logabsdet_jacobian(refbc_f_tpl, xs)
                    @test ly == ly_ref
                    @test ladj ≈ ladj_ref
                end              
            end
        end
    end

    @testset "fchain" begin
        @test @inferred(fchain()) isa FunctionChain
        @test @inferred(fchain(log)) isa FunctionChain
        @test @inferred(fchain(exp)) isa FunctionChain
        @test @inferred(fchain(log, exp)) isa FunctionChain
        @test @inferred(fchain((log, exp))) isa FunctionChain
        @test @inferred(fchain(a = log, c = exp)) isa FunctionChain
        @test @inferred(fchain((a = log, c = exp))) isa FunctionChain
        @test @inferred(fchain([log, exp])) isa FunctionChain
        @test @inferred(fchain(log, ForwardDiff.Dual)) isa FunctionChain
        @test @inferred(fchain(log, AbstractFloat)) isa FunctionChain
        @test @inferred(fchain(ForwardDiff.Dual, fbcast(log))) isa FunctionChain
        @test @inferred(fchain(AbstractFloat, fbcast(log))) isa FunctionChain
        @test @inferred(fchain((Mul(i) for i in 3:4))) isa FunctionChain

        @test_deprecated fchain((log, AbstractFloat)) isa FunctionChain
        @test_throws ArgumentError fchain(a = log, c = AbstractFloat)
        @test_throws ArgumentError fchain((a = log, c = AbstractFloat))

        # Not inferrable:
        @test_deprecated fchain((log, ForwardDiff.Dual)) isa FunctionChain
    end


    @testset "composition" begin
        for f in [identity, fchain(log, exp), fchain((Mul(i) for i in 3:4)), ForwardDiff.Dual, AbstractFloat]
            for g in [identity, fchain(log, exp), fchain((Mul(i) for i in 3:4)), ForwardDiff.Dual, AbstractFloat]
                if f isa FunctionChain || g isa FunctionChain
                    x = 0.73
                    try
                        @test @inferred(f ∘ g) isa FunctionChain
                        @test @inferred((f ∘ g)(x)) == f(g(x))
                    catch
                        global g_state = (;f, g)
                    end
                    @test @inferred(f ∘ g) isa FunctionChain
                    @test @inferred((f ∘ g)(x)) == f(g(x))
                end
            end
        end

        @test @inferred(fchain(Mul(1), Mul(2)) ∘ fchain(Mul(3), Mul(4))) == FunctionChain((Mul(3), Mul(4), Mul(1), Mul(2)))
        @test @inferred(fchain(Mul(3), Mul(4)) ∘ fchain(Mul(1), Mul(2))) == FunctionChain((Mul(1), Mul(2), Mul(3), Mul(4)))
        @test @inferred(fchain([Mul(1), Mul(2)]) ∘ fchain([Mul(3), Mul(4)])) == FunctionChain([Mul(3), Mul(4), Mul(1), Mul(2)])
        @test @inferred(fchain([Mul(3), Mul(4)]) ∘ fchain([Mul(1), Mul(2)])) == FunctionChain([Mul(1), Mul(2), Mul(3), Mul(4)])

        @test @inferred(fchain(Mul(1), Mul(2)) ∘ Mul(3)) == FunctionChain((Mul(3), Mul(1), Mul(2)))
        @test @inferred(Mul(3) ∘ fchain(Mul(1), Mul(2))) == FunctionChain((Mul(1), Mul(2), Mul(3)))
        @test @inferred(fchain([Mul(1), Mul(2)]) ∘ Mul(3)) == FunctionChain([Mul(3), Mul(1), Mul(2)])
        @test @inferred(Mul(3) ∘ fchain([Mul(1), Mul(2)])) == FunctionChain([Mul(1), Mul(2), Mul(3)])

        @test @inferred(fchain(a = sin, x = cos, c = tan) ∘ fchain(e = log, d = exp)) == FunctionChain((e = log, d = exp, a = sin, x = cos, c = tan))
        @test @inferred(fchain(a = sin, x = cos, c = tan) ∘ fchain(e = log, x = exp)) == FunctionChain((FunctionChain((e = log, x = exp)), FunctionChain((a = sin, x = cos, c = tan))))
    end

    let f = fchain(Mul(rand(3,3)), Add(rand(3)), Mul(rand(3,3)))
        @test f == deepcopy(f)

        @static if VERSION >= v"1.9"
            @test @inferred(Adapt.adapt(Array{Float32}, f)) ≈ f
        else
            @test Adapt.adapt(Array{Float32}, f) ≈ f
        end
        @test @inferred(Adapt.adapt(Array{Float32}, f)(rand(Float32, 3))) isa Vector{Float32}

        params, f_ctor = @inferred Functors.functor(f)
        @test @inferred(f_ctor(params)) == f
        @test Functors.fmap(Array{Float32}, f) ≈ f
        @test @inferred(Functors.fmap(Array{Float32}, f)(rand(Float32, 3))) isa Vector{Float32}
    end

    x = 0.3
    xs = [0.3, 0.4, 0.5]

    vx = [1.1, 2.2]
    vxs = [1.1, 2.2], [3.3, 4.4], [5.5, 6.6]

    @test @inferred(fchain(Mul(3), Add(2), InvMul(2.5), Subtract(4))) isa FunctionChain
    f = fchain(Mul(3), Add(2), InvMul(2.5), Subtract(4))
    @test @inferred(f([1.1, 2.2])) == 2.5 \ (3 * [1.1, 2.2] .+ 2) .- 4
    test_function_chain(f, vx, vxs, true, true, "Tuple of AffineMap")
    @test @inferred(((fc)->(fc...,))(f)) == fchainfs(f)

    @test @inferred(fchain(a = Mul(3), c = Add(2), b = InvMul(2.5), z = Subtract(4))) isa FunctionChain
    f = fchain(a = Mul(3), c = Add(2), b = InvMul(2.5), z = Subtract(4))
    @test @inferred(f([1.1, 2.2])) == 2.5 \ (3 * [1.1, 2.2] .+ 2) .- 4
    test_function_chain(f, vx, vxs, true, true, "NamedTuple of AffineMap")
    @test @inferred(merge((;), f)) == fchainfs(f)
    @test (; f...) == fchainfs(f)

    @test @inferred(fchain(Mul.(2:5))) isa FunctionChain{Vector{Mul{Int}}}
    f = fchain(Mul.(2:5))
    @test @inferred(f(1.2)) ≈ 1.2 * prod(2:5)
    test_function_chain(f, x, xs, true, true, "Vector of Mul, scalar arg")
    @test @inferred(with_intermediate_results(f, [1.1, 2.2])) ≈ [[2.2, 4.4], [6.6, 13.2], [26.4, 52.8], [132.0, 264.0]]
    test_function_chain(f, vx, vxs, true, true, "Tuple of Mul, vector arg")

    @test @inferred(fchain(Mul(i) for i in 2:5)) isa FunctionChain{<:Base.Generator}
    f = fchain(Mul(i) for i in 2:5)
    @test @inferred(with_intermediate_results(f, [1.1, 2.2])) ≈ [[2.2, 4.4], [6.6, 13.2], [26.4, 52.8], [132.0, 264.0]]
    test_function_chain(f, vx, vxs, false, true, "Generator of Mul")

    @test @inferred(fchain(repeated(Mul(2), 5))) isa FunctionChain{<:Base.Iterators.Take}
    f = fchain(repeated(Mul(2), 5))
    test_function_chain(f, vx, vxs, true, true, "Generator of Mul")

    @test @inferred(inv ∘ fchain() ∘ (expm1 ∘ log)) == fchain((log, expm1, inv))
    @test @inferred(fchain((log, expm1, inv))(0.3)) == (inv ∘ expm1 ∘ log)(0.3)
    test_function_chain(fchain((log, expm1, inv)), x, xs, true, true, "log-expm1-inv")
    test_function_chain(fchain((a = log, c = expm1, d = inv)), x, xs, true, true, "log-expm1-inv")
    test_function_chain(fchain(fill(expm1, 3)), x, xs, true, true, "fill-expm1")
    test_function_chain(fchain((log, sin, expm1, inv)), x, xs, false, false, "log-sin-expm1-inv")
    test_function_chain(fchain(fill(sin, 3)), x, xs, false, false, "fill-sin")

    @static if isdefined(Main, :FlexiMaps)
        @testset "FlexiMaps support" begin
            @test @inferred(FlexiMaps.islinear(fchain(Mul(4), Mul(3)))) == true
            @test @inferred(FlexiMaps.islinear(fchain(Mul(4), Add(3)))) == false
            @test @inferred(FlexiMaps.islinear(fchain(sin, cos))) == false

            @test @inferred(FlexiMaps.isaffine(fchain(Mul(4), Mul(3)))) == true
            @test @inferred(FlexiMaps.isaffine(fchain(Mul(4), Add(3)))) == true
            @test @inferred(FlexiMaps.isaffine(fchain(sin, cos))) == false
        end
    end

    @testset "ffchain" begin
        @test @inferred(ffchain()) === identity
        @test @inferred(ffchain(identity)) === identity
        @test @inferred(ffchain(identity, identity)) === identity
        @test @inferred(ffchain(log)) === log
        @test @inferred(ffchain(log, exp)) === fchain(log, exp)
        @test @inferred(ffchain(log, exp, sqrt)) === fchain(log, exp, sqrt)
        @test @inferred(ffchain(Int)) isa FunctionChains.AsFunction{Type{Int}}
        @test @inferred(ffchain(identity, Int, identity)) isa FunctionChains.AsFunction{Type{Int}}
        @test @inferred(ffchain(Int, Float32)) === fchain(Int, Float32)
        @test @inferred(ffchain(identity, Float32 ∘ identity ∘ Int, identity)) === fchain(Int, Float32)
        @test @inferred(ffchain((identity ∘ identity) ∘ identity ∘ (identity ∘ identity))) === identity
        @test @inferred(ffchain((sin ∘ cos) ∘ identity ∘ (tan ∘ identity))) === fchain(tan, cos, sin)
        @test @inferred(ffchain((sin ∘ cos) ∘ identity ∘ (tan ∘ identity), fchain(exp, log, sqrt), Float32)) === fchain(tan, cos, sin, exp, log, sqrt, Float32)
    end
end

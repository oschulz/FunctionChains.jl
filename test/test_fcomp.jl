# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test

using AffineMaps
import ForwardDiff, Zygote

@testset "fcomp" begin
    fs_tpl = (log, sqrt, exp)
    fs_array = [Mul(3), Mul(4), Mul(2)]

    x = 0.73

    @test @inferred(fcomp(fs_tpl)) == fchain(reverse(fs_tpl))
    @test @inferred(fcomp(fs_tpl...)) == fchain(reverse(fs_tpl))
    @test @inferred(fcomp(fs_array)) == fchain(reverse(fs_array))
    
    @test @inferred(fcomp(fs_tpl)(x)) == (log ∘ sqrt ∘ exp)(x)
    @test @inferred(fcomp(fs_tpl...)(x)) == (log ∘ sqrt ∘ exp)(x)

    @test_throws ArgumentError fcomp((String, string))

    @testset "ffcomp" begin
        @test @inferred(ffcomp()) === ffchain()
        @test @inferred(ffcomp(identity)) === ffchain(identity)
        @test @inferred(ffcomp(identity, identity)) === ffchain(identity, identity)
        @test @inferred(ffcomp(log)) === ffchain(log)
        @test @inferred(ffcomp(log, exp)) === ffchain(exp, log)
        @test @inferred(ffcomp(log, exp, sqrt)) === ffchain(sqrt, exp, log)
        @test @inferred(ffcomp(Int)) === ffchain(Int)
        @test @inferred(ffcomp(identity, Int, identity)) === ffchain(identity, Int, identity)
        @test @inferred(ffcomp(Int, Float32)) === ffchain(Float32, Int)
        @test @inferred(ffcomp((identity ∘ identity) ∘ identity ∘ (identity ∘ identity))) === identity
        @test @inferred(ffcomp(identity, Float32 ∘ identity ∘ Int, identity)) === ffchain(identity, Float32 ∘ identity ∘ Int, identity)
        @test @inferred(ffcomp((sin ∘ cos) ∘ identity ∘ (tan ∘ identity))) === ffchain((sin ∘ cos) ∘ identity ∘ (tan ∘ identity))
        @test @inferred(ffcomp((sin ∘ cos) ∘ identity ∘ (tan ∘ identity), fchain(exp, log, sqrt), Float32)) === ffchain(Float32, fchain(exp, log, sqrt), (sin ∘ cos) ∘ identity ∘ (tan ∘ identity))

        @test @inferred(ffcomp(FCTestScale(3.0), FCTestScale(2.0))) === FCTestScale(6.0)
        @test @inferred(ffcomp(log10 ∘ exp10)) === identity
        @test @inferred(ffcomp(fchain(FCTestScale(2.0), FCTestScale(3.0)))) === FCTestScale(6.0)
        @test @inferred(ffcomp(fchain(FCTestScale(2.0), fchain(FCTestScale(3.0), FCTestScale(4.0))))) === FCTestScale(24.0)
        @test @inferred(ffcomp(exp, log10 ∘ exp10)) === exp
    end

    @testset "AD through composition construction" begin
        # Zygote must handle compositions constructed inside the
        # differentiated function, including parameters captured in
        # composed closures:
        for compose in (fcomp, ffcomp)
            f_grad2 = x -> compose(y -> y + x, y -> 2 * y)(1.0)
            @test Zygote.gradient(f_grad2, 3.0)[1] == 1.0
            f_grad3 = x -> compose(y -> y + x, y -> 2 * y, y -> y - 1)(1.0)
            @test Zygote.gradient(f_grad3, 3.0)[1] == 1.0
            f_nograd = x -> compose(sin, cos)(x)
            @test Zygote.gradient(f_nograd, 0.3)[1] ≈ ForwardDiff.derivative(f_nograd, 0.3)
        end
    end

end

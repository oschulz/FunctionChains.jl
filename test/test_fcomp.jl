# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test

using AffineMaps

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
    end
end

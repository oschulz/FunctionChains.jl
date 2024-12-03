# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test

using AffineMaps

@testset "fprod" begin
    fs_tpl = (log, sqrt, exp)
    fs_array = [Mul(3), Mul(4), Mul(2)]

    x = 0.73

    @test @inferred(fprod(fs_tpl)) == fchain(reverse(fs_tpl))
    @test @inferred(fprod(fs_tpl...)) == fchain(reverse(fs_tpl))
    @test @inferred(fprod(fs_array)) == fchain(reverse(fs_array))
    
    @test @inferred(fprod(fs_tpl)(x)) == (log ∘ sqrt ∘ exp)(x)
    @test @inferred(fprod(fs_tpl...)(x)) == (log ∘ sqrt ∘ exp)(x)

    @test_throws ArgumentError fprod((String, string))
end

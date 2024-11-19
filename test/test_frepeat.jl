# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test

using InverseFunctions: inverse

@testset "frepeat" begin
    f = sqrt
    n = 5
    fc = frepeat(f, n)
    x = 0.73
    y_ref = f(f(f(f(f(x)))))

    @test @inferred(fc(x)) == y_ref
    @test @inferred(inverse(fc)(y_ref)) ≈ x

    @test @inferred(inverse(inverse(fc))) == fc

    @test f∘̂ 0 == identity
    @test f∘̂ n == fc
    @test f∘̂ -n == inverse(fc)
end

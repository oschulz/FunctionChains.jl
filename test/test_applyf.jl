# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test


@testset "applyf" begin
    foo() = 0
    foo(args...) = sum(args)

    @test @inferred(applyf(foo)) == foo()
    @test @inferred(applyf(foo, 2)) == foo(2)
    @test @inferred(applyf(foo, 2, 4)) == foo(2, 4)
    @test @inferred(applyf(foo, 2, 4, 3)) == foo(2, 4, 3)
end

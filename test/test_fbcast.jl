# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test


@testset "fbcast" begin
    @test @inferred(fbcast(log)) === Base.Broadcast.BroadcastFunction(log)
    @test @inferred(fbcast(Real)) === Base.Broadcast.BroadcastFunction{Type{Real}}(Real)
end

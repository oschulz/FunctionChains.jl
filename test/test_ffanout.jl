# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test

using AffineMaps

@testset "ffanout" begin
    function test_function_fanout(ff, x, label::AbstractString)
        @testset "$label" begin
            @test @inferred(ffanoutfs(ff)) === getfield(ff, :_fs)
            
            @test (ff...,) == (ffanoutfs(ff)...,)
            @test [ff...] == [ffanoutfs(ff)...]

            rf = x -> map(f -> f(x), ffanoutfs(ff))

            @test @inferred(ff(x)) == rf(x)
        end
    end

    fs_tpl = (log, log, exp)
    fs_nt = (a = log, z = log, b = exp)
    fs_vector = [Mul(3), Mul(4), Mul(2)]
    fs_array = Mul.(rand(2,3,4))

    fs_nested = (
        tpl = ffanout(fs_tpl),
        nt = ffanout(fs_nt),
        v = ffanout(fs_vector),
        a = ffanout(fs_array),
    )

    x = rand()

    test_function_fanout(ffanout(fs_tpl), x, "Tuple of functions")
    test_function_fanout(ffanout(fs_tpl...), x, "Splat tuple of functions")
    test_function_fanout(ffanout(fs_nt), x, "NamedTuple of functions")
    test_function_fanout(ffanout(fs_vector), x, "Vector of functions")
    test_function_fanout(ffanout(fs_array), x, "Array of functions")
    test_function_fanout(ffanout(fs_nested), x, "Array of functions")
end

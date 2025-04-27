# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using FunctionChains
using Test

using InverseFunctions: inverse, NoInverse, test_inverse
using ChangesOfVariables: with_logabsdet_jacobian, NoLogAbsDetJacobian, test_with_logabsdet_jacobian
using Accessors: set, PropertyLens, IndexLens

using AffineMaps

include("approx_cmp.jl")
include("getjacobian.jl")
include("modify_output.jl")

@testset "fcprod" begin
    function test_function_product(fp, x, invertible::Bool, settable::Bool, with_ladj::Bool, label::AbstractString)
        @testset "$label" begin
            @test @inferred(fcprodfs(fp)) === getfield(fp, :_fs)
            
            @test (fp...,) == (fcprodfs(fp)...,)
            @test [fp...] == [fcprodfs(fp)...]

            rf = x -> map(applyf, fcprodfs(fp), x)

            @test @inferred(fp(x)) == rf(x)
            y_out = fp(x)

            if invertible
                @test @inferred(inverse(fp)) isa FCartProd
                test_inverse(fp, x, compare = approx_cmp)
            else
                @test @inferred(inverse(fp)) isa NoInverse{typeof(fp)}
            end

            if settable
                new_y = modify_output(y_out)
                @test @inferred(set(x, fp, new_y)) == map(set, x, fp._fs, new_y)
            end

            if with_ladj
                y_ladj_ref = rf(x), sum(map(x -> x[2], map(with_logabsdet_jacobian, fcprodfs(fp), x)))
                @test approx_cmp(@inferred(with_logabsdet_jacobian(fp, x)), y_ladj_ref)
                test_with_logabsdet_jacobian(fp, x, getjacobian; compare = approx_cmp)
            else
                @test @inferred(with_logabsdet_jacobian(fp, x)) isa NoLogAbsDetJacobian{typeof(fp), typeof(x)}
            end
        end
    end


    fs_tpl = (log, log, exp)
    fs_nt = (a = log, z = log, b = exp)
    fs_vector = [Mul(3), Mul(4), Mul(2)]
    fs_array = Mul.(rand(2,3,4))

    fs_nested = (
        tpl = fcprod(fs_tpl),
        nt = fcprod(fs_nt),
        v = fcprod(fs_vector),
        a = fcprod(fs_array),
    )

    fs_tpl_noinv = (log, sin, exp)
    fs_nt_noinv = (a = log, z = sin, b = exp)
    fs_vector_noinv = fill(sin, 3)
    fs_array_noinv = fill(sin, 2,3,4)

    fs_nested_noinv = (
        tpl = fcprod(fs_tpl_noinv),
        nt = fcprod(fs_nt_noinv),
        v = fcprod(fs_vector_noinv),
        a = fcprod(fs_array_noinv),
    )

    fs_tpl_identity = (identity, identity, identity)
    fs_nt_identity = (a = identity, z = identity, b = identity)
    fs_vector_identity = fill(identity, 3)
    fs_array_identity = fill(identity, 2,3,4)

    fs_nested_identity = (
        tpl = fcprod(fs_tpl_identity),
        nt = fcprod(fs_nt_identity),
        v = fcprod(fs_vector_identity),
        a = fcprod(fs_array_identity),
    )

    x_tpl = (rand(), rand(), rand())
    x_nt = (a = rand(), z = rand(), b = rand())
    x_vector = [rand(), rand(), rand()]
    x_array = rand(2,3,4)

    x_nested = (
        tpl = x_tpl,
        nt = x_nt,
        v = x_vector,
        a = x_array,    
    )


    test_function_product(fcprod(fs_tpl), x_tpl, true, true, true, "Tuple of functions")
    test_function_product(fcprod(fs_tpl...), x_tpl, true, true, true, "Splat tuple of functions")
    test_function_product(fcprod(fs_nt), x_nt, true, true, true, "NamedTuple of functions")
    test_function_product(fcprod(fs_vector), x_vector, true, true, true, "Vector of functions")
    test_function_product(fcprod(fs_array), x_array, true, true, true, "Array of functions")
    test_function_product(fcprod(fs_nested), x_nested, true, true, true, "Array of functions")

    test_function_product(fcprod(fs_tpl_noinv), x_tpl, false, false, false, "Non-Inv Tuple of functions")
    test_function_product(fcprod(fs_tpl_noinv...), x_tpl, false, false, false, "Non-Inv Splat tuple of functions")
    test_function_product(fcprod(fs_nt_noinv), x_nt, false, false, false, "Non-Inv NamedTuple of functions")
    test_function_product(fcprod(fs_vector_noinv), x_vector, false, false, false, "Non-Inv Vector of functions")
    test_function_product(fcprod(fs_array_noinv), x_array, false, false, false, "Non-Inv Array of functions")
    test_function_product(fcprod(fs_nested_noinv), x_nested, false, false, false, "Non-Inv Array of functions")

    test_function_product(fcprod(fs_tpl_identity), x_tpl, true, true, true, "Non-Inv Tuple of functions")
    test_function_product(fcprod(fs_tpl_identity...), x_tpl, true, true, true, "Non-Inv Splat tuple of functions")
    test_function_product(fcprod(fs_nt_identity), x_nt, true, true, true, "Non-Inv NamedTuple of functions")
    test_function_product(fcprod(fs_vector_identity), x_vector, true, true, true, "Non-Inv Vector of functions")
    test_function_product(fcprod(fs_array_identity), x_array, true, true, true, "Non-Inv Array of functions")
    test_function_product(fcprod(fs_nested_identity), x_nested, true, true, true, "Non-Inv Array of functions")

    test_function_product(
        fcprod(
            fchain(PropertyLens{:d}(), exp),
            sqrt,
            IndexLens(2)
        ),
        (
            (d = 4.2, c = "Hello"),
            42,
            [3, 5, 2]
        ),
        false, true, false, "fcprod with lenses"
    )

    test_function_product(
        fcprod((
            c = fchain(PropertyLens{:d}(), exp),
            d = sqrt,
            e = IndexLens(2)
        )),
        (
            c = (d = 4.2, c = "Hello"),
            d = 42,
            e = [3, 5, 2]
        ),
        false, true, false, "fcprod with lenses"
    )

    test_function_product(
        fcprod([log ∘ IndexLens(i) for i in[2,4,5]]), 
        [rand(10) for i in 1:3],
        false, true, false, "fcprod with lenses"
    )

    let f = fcprod(Mul(rand(3,3)), Add(rand(3)), Mul(rand(3,3)))
        @test f == deepcopy(f)
        @test f ≈ deepcopy(f)

        x = (rand(Float32, 3),rand(Float32, 3),rand(Float32, 3))
        @test @inferred(Adapt.adapt(Array{Float32}, f)) ≈ f
        @test @inferred(Adapt.adapt(Array{Float32}, f)(x)) isa NTuple{3,Vector{Float32}}

        params, f_ctor = @inferred Functors.functor(f)
        @test @inferred(f_ctor(params)) == f
        @test Functors.fmap(Array{Float32}, f) ≈ f
        @test @inferred(Functors.fmap(Array{Float32}, f)(x)) isa NTuple{3,Vector{Float32}}
    end
end

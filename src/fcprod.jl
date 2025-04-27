# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    fcprodfs(fp)

Get the component functions of a cartesian product of functions.

See [`fcprod`](@ref) for details.
"""
function fcprodfs end
export fcprodfs


"""
    struct FCartProd{FS}<:Function

Represents a
[Cartesian product of functions](https://en.wikipedia.org/wiki/Cartesian_product#Cartesian_product_of_functions).

Use``
A `FCartProd` has a single field `fs` which may be a `Tuple`, `NamedTuple`, an
array or a generator/iterator of functions.

Use [`fcprod`](@ref) to construct products of functions instead of using the
constructor `FCartProd(fs)` directly.
"""
struct FCartProd{FS} <: Function
    _fs::FS
end
export FCartProd

fcprodfs(fp::FCartProd) = getfield(fp, :_fs)

function Base.show(io::IO, m::MIME"text/plain", fp::FCartProd{<:Tuple})
    print(io, "fcprod")
    show(io, m, fp._fs)
end

Base.show(io::IO, fp::FCartProd) = show(io, MIME"text/plain"(), fp)

# Enables `(; fp...)`
Base.merge(a::NamedTuple, fp::FCartProd{<:NamedTuple{names}}) where {names} = merge(a, fp._fs)

Base.merge(a::FCartProd{<:NamedTuple}) = a
function Base.merge(a::FCartProd{<:NamedTuple}{names}, b::FCartProd{<:NamedTuple}, cs::FCartProd{<:NamedTuple}...) where {names}
    return merge(fcprod(; a..., b...), cs...)
end

# ToDo: Define length?
# Base.length(fp::FCartProd) = length(fp._fs)

# Enables `(fp...,)` and `[fp...]`
Base.iterate(fp::FCartProd) = iterate(fp._fs)
Base.iterate(fp::FCartProd, state) = iterate(fp._fs, state)
    

(fp::FCartProd{<:Tuple{Vararg{Any,N}}})(x::Tuple{Vararg{Any,N}}) where N = map(applyf, fp._fs, x)

(::FCartProd{<:Tuple{Vararg{typeof(identity),N}}})(x::Tuple{Vararg{Any,N}}) where N = x

function (@nospecialize(fs::FCartProd{<:Tuple}))(@nospecialize(x::Tuple))
    throw(ArgumentError("Can't apply FCartProd over Tuple of length $(length(fp._fs)) to Tuple of length $(length(x))."))
end


function (fp::FCartProd{<:NamedTuple{names}})(x::NamedTuple{names}) where names
    NamedTuple{names}(map(applyf, values(fp._fs), values(x)))
end

(::FCartProd{<:NamedTuple{names,<:Tuple{Vararg{typeof(identity)}}}})(x::NamedTuple{names}) where names = x

function (@nospecialize(fs::FCartProd{<:NamedTuple}))(@nospecialize(x::NamedTuple))
    throw(ArgumentError("Can't apply FCartProd over NamedTuple with names $(propertynames(fp._fs)) to NamedTuple with names $(propertynames(x))."))
end


Base.@propagate_inbounds function (fp::FCartProd)(x)
    fs = fp._fs
    @boundscheck _check_fp_sizes(IteratorSize(fs), IteratorSize(x), fs, x)
    @inbounds result = _fp_apply(fs, x)
    return result
end

_check_fp_sizes(::IteratorSize, ::IteratorSize, ::Any, ::Any) = nothing

function _check_fp_sizes(::HasShape, ::HasShape, fs, x)
    size_fs, size_x = size(fs), size(x)
    if size_fs != size_x
        throw(ArgumentError("Can't apply FCartProd of size $size_fs to argument of size $size_x."))
    end
    return nothing
end

function _check_fp_sizes(::HasLength, ::HasLength, fs, x)
    length_fs, length_x = length(fs), length(x)
    if length_fs != length_x
        throw(ArgumentError("Can't apply FCartProd of length $length_fs to argument of length $length_x."))
    end
    return nothing
end

Base.@propagate_inbounds _fp_apply(fs, x) = applyf.(fs, x)

_fp_apply(::AbstractArray{typeof(identity)}, x) = x


function InverseFunctions.inverse(fp::FCartProd)
    inv_fs = map(inverse, fp._fs)
    if _contains_noinverse(inv_fs) isa Val{true}
        NoInverse(fp)
    else
        FCartProd(inv_fs)
    end
end

# Mapping directly over fp._fs::NamedTuple causes type instability with
# NoInverse for some reason, so need specialized method for FCartProd{<:NamedTuple}:
function InverseFunctions.inverse(fp::FCartProd{<:NamedTuple{names}}) where names
    inv_fs = map(inverse, values(fp._fs))
    if _contains_noinverse(inv_fs) isa Val{true}
        NoInverse(fp)
    else
        FCartProd(NamedTuple{names}(inv_fs))
    end
end


"""
    fcprod()
    fcprod(fs)
    fcprod(fs...)

Construct the
[Cartesian product of functions](https://en.wikipedia.org/wiki/Cartesian_product#Cartesian_product_of_functions)
over the functions `fs`.

Typically returns a [`FCartProd`](@ref), but may be specialized to
return other product function types for specific types of element functions.

`fs` must be iterable, it may be a `Tuple`, `NamedTuple`, an array or a
generator/iterator of functions.

`fcprod` behaves like

```julia
fcprod((f_a, f_b, ...))((x_a, x_b, ...)) = (f_a(x_a), f_b(x_b), ...)
fcprod((a = f_a, b = f_b, ...))((a = x_a, b = x_b, ...)) = (a = f_a(x_a), b = f_b(x_b), ...)
fcprod([f_a, f_b, ...])([x_a, x_b, ...]) = [f_a(x_a), f_b(x_b, ...)]
```

This is similar, semantically, to Haskell's `***` for arrows.

For `fp = fcprod(fs)`, use [`fcprodfs(fp)`](@ref) to retrieve `fs`.

The resulting product of functions supports `InverseFunctions.inverse` and/or
`ChangesOfVariables.with_logabsdet_jacobian` if all functions in the product
do so.
"""
function fcprod end
export fcprod

@inline fcprod(;fs...) = _fprod_vararg_impl(values(fs))
@inline _fprod_vararg_impl(::NamedTuple{()}) = fcprod(())
@inline _fprod_vararg_impl(fs::NamedTuple) = fcprod(fs)

@inline fcprod(fs::FS) where FS = _fprod_onearg(fs, Val(static_hasmethod(iterate, Tuple{FS})))
_fprod_onearg(fs::FS, ::Val{true}) where FS = FCartProd(fs)
_fprod_onearg(f::F, ::Val{false}) where F = FCartProd(_typed_funcs_tuple(f))

@inline fcprod(fs::Vararg{Any}) = FCartProd(_typed_funcs_tuple(fs...))


@inline fcprod(fs::Tuple{Vararg{Function}}) = FCartProd(fs)

@inline function fcprod(fs::Tuple)
    _check_fp_contents(fs)
    FCartProd(fs)
end

function _check_fp_contents(fs::Tuple)
    if any(Base.Fix2(isa, Type), fs)
        throw(ArgumentError("Do not use fcprod(fs::Tuple) with fs elements that are types, due to possible type instabilities, use `fcprod(fs...)` instead."))
    end
end


@inline fcprod(fs::NamedTuple{names,<:Tuple{Vararg{Function}}}) where names = FCartProd(fs)

@inline function fcprod(fs::NamedTuple)
    _check_fp_contents(fs)
    FCartProd(fs)
end

function _check_fp_contents(fs::NamedTuple)
    if any(Base.Fix2(isa, Type), values(fs))
        throw(ArgumentError("Do not use fcprod(fs::NamedTuple) with fs elements that are types, due to possible type instabilities."))
    end
end

# This file is a part of FunctionProducts.jl, licensed under the MIT License (MIT).


#=
struct FunctionProduct{FS<:AbstractVector{<:Function}} <: Function
    fs::FS
end

(f::FunctionProduct)(x::AbstractVector) = x .|> f.fs

InverseFunctions.inverse(f::FunctionProduct) = FunctionProduct(inverse.(f.fs))

function ChangesOfVariables.with_logabsdet_jacobian(f::FunctionProduct, x::AbstractVector)
    ys_ladjs = with_logabsdet_jacobian.(f.fs, x)
    ys = map(x -> x[1], ys_ladjs)
    ladjs = map(x -> x[2], ys_ladjs)
    return ys, sum(ladjs)
end

=#


"""
    struct FunctionProduct{FS}<:Function

Represents Cartesian product of functions.

A `FunctionProduct` has a single field `fs` which may be a `Tuple`,
`NamedTuple`, an array or a generator/iterator of functions.

Use [`fprod`](@ref) to construct products of functions instead of using the
constructor `FunctionProduct(fs)` directly.
"""
struct FunctionProduct{FS} <: Function
    fs::FS
end
export FunctionProduct

function Base.show(io::IO, m::MIME"text/plain", fp::FunctionProduct{<:Tuple})
    print(io, "fprod")
    show(io, m, fp.fs)
end

Base.show(io::IO, fp::FunctionProduct) = show(io, MIME"text/plain"(), fp)

# Enables `(; fp...)`
Base.merge(a::NamedTuple, fp::FunctionProduct{<:NamedTuple{names}}) where {names} = merge(a, fp.fs)

Base.merge(a::FunctionProduct{<:NamedTuple}) = a
Base.merge(a::FunctionProduct{<:NamedTuple}{names}, b::FunctionProduct{<:NamedTuple}, cs::FunctionProduct{<:NamedTuple}...) where {names} =
    merge(fprod(; a..., b...), cs...)


(fp::FunctionProduct{<:Tuple{Vararg{Any,N}}})(x::Tuple{Vararg{Any,N}}) where N = map(applyf, fp.fs, x)

(::FunctionProduct{<:Tuple{Vararg{typeof(identity),N}}})(x::Tuple{Vararg{Any,N}}) where N = x

function (@nospecialize(fs::FunctionProduct{<:Tuple}))(@nospecialize(x::Tuple))
    throw(ArgumentError("Can't apply FunctionProduct over Tuple of length $(length(fp.fs)) to Tuple of length $(length(x))."))
end


function (fp::FunctionProduct{<:NamedTuple{names}})(x::NamedTuple{names}) where names
    NamedTuple{names}(map(applyf, values(fp.fs), values(x)))
end

(::FunctionProduct{<:NamedTuple{names,<:Tuple{Vararg{typeof(identity)}}}})(x::NamedTuple{names}) where names = x

function (@nospecialize(fs::FunctionProduct{<:NamedTuple}))(@nospecialize(x::NamedTuple))
    throw(ArgumentError("Can't apply FunctionProduct over NamedTuple with names $(propertynames(fp.fs)) to NamedTuple with names $(propertynames(x))."))
end


Base.@propagate_inbounds function (fp::FunctionProduct)(x)
    @boundscheck _check_fp_sizes(fp, x)
    @inbounds result = _fp_apply(fp, x)
    return result
end

function _check_fp_sizes(fp::FunctionProduct, x)
    if size(fp.fs) != size(x)
        throw(ArgumentError("Can't apply FunctionProduct of size $(size(fp.fs)) to argument of size $(size(x))."))
    end
    return nothing
end

Base.@propagate_inbounds _fp_apply(fp::FunctionProduct, x) = applyf.(fp.fs, x)

_fp_apply(::FunctionProduct{<:AbstractArray{typeof(identity)}}, x) = x


"""
    fprod()
    fprod(fs)
    fprod(fs...)

Construct a Cartesian product over the functions `fs`.

Typically returns a [`FunctionProduct`](@ref), but may be specialized to
return other product function types for specific types of element functions.

`fs` must be iterable, it may be a `Tuple`, `NamedTuple`, an array or a
generator/iterator of functions.

`fprod` behaves like

```julia
fprod((f_a, f_b, ...))((x_a, x_b, ...)) = (f_a(x_a), f_b(x_b), ...)
fprod((a = f_a, b = f_b, ...))((a = x_a, b = x_b, ...)) = (a = f_a(x_a), b = f_b(x_b), ...)
fprod([f_a, f_b, ...])([x_a, x_b, ...]) = [f_a(x_a), f_b(x_b, ...)]
```

The resulting product of functions supports `InverseFunctions.inverse` and/or
`ChangesOfVariables.with_logabsdet_jacobian` if all functions in the product
do so.
"""
function fprod end
export fprod

@inline fprod(;fs...) = _fprod_vararg_impl(values(fs))
@inline _fprod_vararg_impl(::NamedTuple{()}) = fprod(())
@inline _fprod_vararg_impl(fs::NamedTuple) = fprod(fs)

@inline fprod(fs::FS) where FS = _fprod_onearg(fs, Val(static_hasmethod(iterate, Tuple{FS})))
_fprod_onearg(fs::FS, ::Val{true}) where FS = FunctionProduct(fs)
_fprod_onearg(f::F, ::Val{false}) where F = FunctionProduct(_typed_funcs_tuple(f))

@inline fprod(fs::Tuple{Vararg{Function}}) = FunctionProduct(fs)

@inline fprod(fs::Vararg{Any}) = FunctionProduct(_typed_funcs_tuple(fs...))

function fprod(@nospecialize(::Tuple))
    throw(ArgumentError("Do not use fprod(fs::Tuple) with fs elements not of type Function, due to possible type instabilities, use `fprod(fs...)` instead."))
end

@inline fprod(fs::NamedTuple{names,<:Tuple{Vararg{Function}}}) where names = FunctionProduct(fs)

function fprod(@nospecialize(::NamedTuple))
    throw(ArgumentError("Do not use fprod(fs::NamedTuple) or fprod(;fs...) with fs elements not of type Function, due to type instability."))
end

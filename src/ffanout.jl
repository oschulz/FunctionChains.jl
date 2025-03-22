# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    ffanoutfs(fc)

Get the component functions of a function fanout.

See [`ffanout`](@ref) for details.
"""
function ffanoutfs end
export ffanoutfs


"""
    struct FFanout{FS}<:Function

Represents a function fanout.

Use``
A `FFanout` has a single field `fs` which may be a `Tuple`, `NamedTuple`, an
array or a generator/iterator of functions.

Use [`ffanout`](@ref) to construct a function fanout instead of using the
constructor `FFanout(fs)` directly.
"""
struct FFanout{FS} <: Function
    _fs::FS
end
export FFanout

ffanoutfs(ff::FFanout) = getfield(ff, :_fs)

function Base.show(io::IO, m::MIME"text/plain", ff::FFanout{<:Tuple})
    print(io, "ffanout")
    show(io, m, ff._fs)
end

Base.show(io::IO, ff::FFanout) = show(io, MIME"text/plain"(), ff)

# Enables `(; ff...)`
Base.merge(a::NamedTuple, ff::FFanout{<:NamedTuple{names}}) where {names} = merge(a, ff._fs)

Base.merge(a::FFanout{<:NamedTuple}) = a
function Base.merge(a::FFanout{<:NamedTuple}{names}, b::FFanout{<:NamedTuple}, cs::FFanout{<:NamedTuple}...) where {names}
    return merge(ffanout(; a..., b...), cs...)
end

# ToDo: Define length?
# Base.length(fc::FFanout) = length(fc._fs)

# Enables `(fc...,)` and `[fc...]`
Base.iterate(fc::FFanout) = iterate(fc._fs)
Base.iterate(fc::FFanout, state) = iterate(fc._fs, state)
    

(ff::FFanout{<:Tuple})(x) = map(Base.Fix2(applyf, x), ff._fs)

function (ff::FFanout{<:NamedTuple{names}})(x) where names
    NamedTuple{names}(map(Base.Fix2(applyf, x), values(ff._fs)))
end

Base.@propagate_inbounds (ff::FFanout)(x) = Base.Fix2(applyf, x).(ff._fs)


"""
    ffanout()
    ffanout(fs)
    ffanout(fs...)

Construct a function fanout with the functions `fs`.

This is similar, semantically, to Haskell's `&&&` for arrows.

Typically returns a [`FFanout`](@ref), but may be specialized to
return other function fanout types for specific types of element functions.

`fs` must be iterable, it may be a `Tuple`, `NamedTuple`, an array or a
generator/iterator of functions.

The resulting function fanout behaves like

```julia
ffanout((f_a, f_b, ...))(x) = (f_a(x), f_b(x), ...)
ffanout((a = f_a, b = f_b, ...))(x) = (a = f_a(x), b = f_b(x), ...)
ffanout([f_a, f_b, ...])(x) = [f_a(x), f_b(x, ...)]
```

For `ff = ffanout(fs)`, use [`ffanoutfs(fc)`](@ref) to retrieve `fs`.
"""
function ffanout end
export ffanout

@inline ffanout(;fs...) = _ffan_vararg_impl(values(fs))
@inline _ffan_vararg_impl(::NamedTuple{()}) = ffanout(())
@inline _ffan_vararg_impl(fs::NamedTuple) = ffanout(fs)

@inline ffanout(fs::FS) where FS = _ffan_onearg(fs, Val(static_hasmethod(iterate, Tuple{FS})))
_ffan_onearg(fs::FS, ::Val{true}) where FS = FFanout(fs)
_ffan_onearg(f::F, ::Val{false}) where F = FFanout(_typed_funcs_tuple(f))

@inline ffanout(fs::Tuple{Vararg{Function}}) = FFanout(fs)

@inline ffanout(fs::Vararg{Any}) = FFanout(_typed_funcs_tuple(fs...))

function ffanout(@nospecialize(::Tuple))
    throw(ArgumentError("Do not use ffanout(fs::Tuple) with fs elements not of type Function, due to possible type instabilities, use `ffanout(fs...)` instead."))
end

@inline ffanout(fs::NamedTuple{names,<:Tuple{Vararg{Function}}}) where names = FFanout(fs)

function ffanout(@nospecialize(::NamedTuple))
    throw(ArgumentError("Do not use ffanout(fs::NamedTuple) or ffanout(;fs...) with fs elements not of type Function, due to type instability."))
end

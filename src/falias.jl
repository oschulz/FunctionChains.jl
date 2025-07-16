# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

"""
    struct FunctionChains.FAlias{F} <: Function

Represents a function alias that behaves like the original function, except
for broadcasting specialization.

User code should not instantiate `FAlias` directly, use [`falias`](@ref)
instead.
"""
struct FAlias{F} <: Function
    _f_wrapped::F
end
# Ensure type stability if function is a type (constructor):
FAlias(::Type{T}) where T = FAlias{Type{T}}(T) # For type stability if 

Base.:(==)(a::FAlias, b::FAlias) = a._f_wrapped == b._f_wrapped
Base.isapprox(a::FAlias, b::FAlias; kwargs...) = isapprox(a._f_wrapped, b._f_wrapped; kwargs...)

@inline (ff::FAlias)(x) = ff._f_wrapped(x)
@inline (ff::FAlias)(xs...; kwargs...) = ff._f_wrapped(xs...; kwargs...)

InverseFunctions.inverse(ff::FAlias) = falias(InverseFunctions.inverse(ff._f_wrapped))


"""
    falias(_f_wrapped)::Function

Wraps a function into an alias that behaves like the original function `_f_wrapped`.

Use [`forig(_f_wrapped)`](@ref) to retrieve the original function:

```julia
g = falias(_f_wrapped)
g(xs...; kwargs...) == _f_wrapped(xs...; kwargs...)

forig(g) == _f_wrapped
```

The resulting function inherits support for Accessors, Adapt,
InverseFunctions, etc., but doesn't inherit broadcasting specialization. This
makes it useful for packages that provide default broadcasting specializations
for certain types of functions.

Example:

See [`ffanout`](@ref) for details.

# Implementation

`falias` will typically return an instance of [`FAlias`](@ref).
"""
function falias end
export falias

@inline falias(_f_wrapped::Function) = _f_wrapped
@inline falias(_f_wrapped) = FAlias(_f_wrapped)
# Ensure type stability if function is a type (constructor):
@inline falias(_f_wrapped::Type{T}) where T = FAlias{Type{T}}(_f_wrapped)


"""
    forig(_f_wrapped)

If `_f_wrapped` is a [`FAlias`](@ref), returns the original function it wraps,
otherwise returns `_f_wrapped`.
"""
function typed_callable end
export typed_callable

@inline typed_callable(_f_wrapped) = _f_wrapped
# Ensure type stability if function is a type (constructor):
@inline typed_callable(_f_wrapped::Type{T}) where T = FAlias{Type{T}}(_f_wrapped)

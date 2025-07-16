# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

"""
    struct FunctionChains.FNoBCSpec{F} <: Function

Represents a function alias that behaves like the original function, except
for broadcasting specialization (in contrast to [`AsFunction`](@ref)).

User code should typically not instantiate `FNoBCSpec` objects directly, but
use [`fnobcspec(f)`](@ref) instead.
"""
struct FNoBCSpec{F} <: Function
    _f_wrapped::F
end
# Ensure type stability if function is a type (constructor):
FNoBCSpec(::Type{T}) where T = FNoBCSpec{Type{T}}(T)

Base.:(==)(a::FNoBCSpec, b::FNoBCSpec) = a._f_wrapped == b._f_wrapped
Base.isapprox(a::FNoBCSpec, b::FNoBCSpec; kwargs...) = isapprox(a._f_wrapped, b._f_wrapped; kwargs...)

@inline (ff::FNoBCSpec)(x) = ff._f_wrapped(x)
@inline (ff::FNoBCSpec)(xs...; kwargs...) = ff._f_wrapped(xs...; kwargs...)

InverseFunctions.inverse(ff::FNoBCSpec) = fnobcspec(InverseFunctions.inverse(ff._f_wrapped))


"""
    fnobcspec(_f_wrapped)::Function

Strips broadcasting specialization from a function.

`fnobcspec` wraps a function into an alias (typically an instance of
[`FunctionChains.FNoBCSpec`](ref)) that behaves like the original function
`_f_wrapped`, except for broadcasting specialization. This makes it
particularly useful for applications that want to provide custom broadcasting
specializations for certain (super-)types of functions that may need to
fall back to standard broadcasting in some cases. Using `falias` in the
fallback path will then prevent infinite recursion.

Use [`forig(_f_wrapped)`](@ref) to retrieve the original function:

```julia
g = fnobcspec(_f_wrapped)
g(xs...; kwargs...) == _f_wrapped(xs...; kwargs...)

forig(g) == _f_wrapped
```

`g` inherits support for Accessors, Adapt, InverseFunctions, etc., from `f`.,but doesn't inherit broadcasting specialization. 
```
"""
function fnobcspec end
export fnobcspec

@inline fnobcspec(_f_wrapped::Function) = _f_wrapped
@inline fnobcspec(_f_wrapped) = FNoBCSpec(_f_wrapped)
# Ensure type stability if function is a type (constructor):
@inline fnobcspec(_f_wrapped::Type{T}) where T = FNoBCSpec{Type{T}}(_f_wrapped)


"""
    forig(_f_wrapped)

If `_f_wrapped` is a [`FNoBCSpec`](@ref), returns the original function it wraps,
otherwise returns `_f_wrapped`.
"""
function forig end

#!!!!!!!!!!!!

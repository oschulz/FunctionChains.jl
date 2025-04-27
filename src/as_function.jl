# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

"""
    struct FunctionChains.AsFunction{F} <: Function

Wraps a callable object to make it a `Function`.

User code should typically not instantiate `AsFunction` objects directly,
but use [`asfunction(f)`](@ref) instead.
"""
struct AsFunction{F} <: Function
    f::F
end
# Ensure type stability if function is a type (constructor):
AsFunction(::Type{T}) where T = AsFunction{Type{T}}(T) # For type stability if 

Base.:(==)(a::AsFunction, b::AsFunction) = a.f == b.f
Base.isapprox(a::AsFunction, b::AsFunction; kwargs...) = isapprox(a.f, b.f; kwargs...)

@inline (ff::AsFunction)(x) = ff.f(x)
@inline (ff::AsFunction)(xs...; kwargs...) = ff.f(xs...; kwargs...)

@inline Broadcast.broadcasted(ff::AsFunction, x) = Broadcast.broadcasted(ff.f, xx)
@inline Broadcast.broadcasted(ff::AsFunction, xs...) = Broadcast.broadcasted(ff.f, xs...)

InverseFunctions.inverse(ff::AsFunction) = _typestable_func(InverseFunctions.inverse(ff.f))

@inline _typestable_func(f) = f
@inline _typestable_func(f::Type{T}) where T = AsFunction{Type{T}}(f)


"""
    asfunction(f)::Function
    asfunction(f::Function) === f

Wraps a callable `f` to make it a `Function`.

If `f isa Function`, simply returns f. If `f` is A
type (constructor), returns a properly typed function
object.
"""
function asfunction end
export asfunction

@inline asfunction(f::Function) = f
@inline asfunction(f) = AsFunction(f)
# Ensure type stability if function is a type (constructor):
@inline asfunction(f::Type{T}) where T = AsFunction{Type{T}}(f)


"""
    typed_callable(f)

Wraps a callable `f` in a typed function object if necessary,
e.g. if `f` is a type (constructor).
"""
function typed_callable end
export typed_callable

@inline typed_callable(f) = f
# Ensure type stability if function is a type (constructor):
@inline typed_callable(f::Type{T}) where T = AsFunction{Type{T}}(f)

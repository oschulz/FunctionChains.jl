# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    with_intermediate_results(f, x)

Apply multi-step function `f` to `x` and return a collection that contains
the intermediate results and the final result.
"""
function with_intermediate_results end
export with_intermediate_results

with_intermediate_results(f, x) = (f(x),)


struct _AsFunction{F} <: Function
    f::F
end

Base.:(==)(a::_AsFunction, b::_AsFunction) = a.f == b.f
Base.isapprox(a::_AsFunction, b::_AsFunction; kwargs...) = isapprox(a.f, b.f; kwargs...)

@inline (ff::_AsFunction{F})(xs...) where F = ff.f(xs...)

@inline _typed_func(f) = f
@inline _typed_func(f::Type{T}) where T = _AsFunction{Type{T}}(f)::_AsFunction{Type{T}}

@inline @generated function _typed_funcs_tuple(fs::Vararg{Any,N}) where N
    expr = Expr(:tuple)
    for i in 1:N
        if fs[i] <: Type
            push!(expr.args, :(_AsFunction{$(fs[i])}(fs[$i])))
        # Future option (breaking) - flatten FunctionChains:
        # elseif fs[i] <: FunctionChain{<:Tuple}
        #    push!(expr.args, :(fs[$i].fs...))
        else
            # Future option (breaking) - remove identities:
            #if !(fs[i] <: typeof(identity))
                push!(expr.args, :(fs[$i]))
            #end
        end
    end
    return expr
end


"""
    struct FunctionChain{FS}<:Function

Represents a chain of composed functions.

A `FunctionChain` has a single field `fs` which may be a tuple, array or
generator/iterator of functions.

`(chain::FunctionChain)(x)` applies the functions in the cain in order of
iteration over `chain.fs`.
```

Supports [`with_intermediate_results`](@ref). Also supports
`InverseFunctions.inverse` and/or `ChangesOfVariables.with_logabsdet_jacobian`
if all functions in the chain do so.

Use [`fchain`](@ref) to construct function chains instead of `FunctionChain(fs)`.
"""
struct FunctionChain{FS} <: Function
    fs::FS
end
export FunctionChain


Base.:(==)(a::FunctionChain, b::FunctionChain) = a.fs == b.fs
Base.isapprox(a::FunctionChain, b::FunctionChain; kwargs...) = _isapprox(a.fs, b.fs; kwargs...)

_isapprox(a, b; kwargs...) = isapprox(a, b; kwargs...)
_isapprox(a::Tuple{Vararg{Any,N}}, b::Tuple{Vararg{Any,N}}; kwargs...) where N = all(map((a, b) ->_isapprox(a, b; kwargs...), a, b))


function Base.show(io::IO, m::MIME"text/plain", fc::FunctionChain)
    print(io, "fchain(")
    show(io, m, fc.fs)
    print(io, ")")
end

function Base.show(io::IO, m::MIME"text/plain", fc::FunctionChain{<:Tuple})
    print(io, "fchain")
    show(io, m, fc.fs)
end

Base.show(io::IO, fc::FunctionChain) = show(io, MIME"text/plain"(), fc)


convert(::Type{FunctionChain}, f::ComposedFunction) = ComposedFunction(_flatten_composed(f))

@inline _flatten_composed(f::F) where {F} = (_typed_func(f),)
@inline _flatten_composed(f::F) where {F<:ComposedFunction} = _typed_funcs_tuple(_flatten_composed(f.inner)..., _flatten_composed(f.outer)...)

@inline Base.:(∘)(f::FunctionChain, g::FunctionChain) = _compose_fc_fc(f, g)
@inline Base.:(∘)(f::FunctionChain, g::ComposedFunction) = _compose_fc_fc(f, FunctionChain(_flatten_composed(g)))
@inline Base.:(∘)(f::ComposedFunction, g::FunctionChain) = _compose_fc_fc(FunctionChain(_flatten_composed(f)), g)
@inline Base.:(∘)(f::FunctionChain, g) = _compose_fc_sf(f, _typed_func(g))
@inline Base.:(∘)(f, g::FunctionChain) = _compose_sf_fc(_typed_func(f), g)
@inline Base.:(∘)(f::FunctionChain, ::typeof(identity)) = f
@inline Base.:(∘)(::typeof(identity), g::FunctionChain) = g

@inline _compose_fc_fc(f::FunctionChain, g::FunctionChain) = FunctionChain((g, f))
@inline _compose_fc_fc(f::FunctionChain{<:Tuple}, g::FunctionChain) = FunctionChain((g, f.fs...))
@inline _compose_fc_fc(f::FunctionChain, g::FunctionChain{<:Tuple}) = FunctionChain((g.fs..., f))
@inline _compose_fc_fc(f::FunctionChain{<:Tuple}, g::FunctionChain{<:Tuple}) = FunctionChain((g.fs..., f.fs...))
@inline _compose_fc_fc(f::FunctionChain{<:AbstractVector{F}}, g::FunctionChain{<:AbstractVector{F}}) where F = FunctionChain(vcat(g.fs, f.fs))

@inline _compose_fc_sf(f::FunctionChain, g) = FunctionChain((g, f))
@inline _compose_fc_sf(f::FunctionChain{<:Tuple}, g) = FunctionChain((g, f.fs...))
@inline _compose_fc_sf(f::FunctionChain{<:AbstractVector{F}}, g::F) where F = FunctionChain(pushfirst!(copy(f.fs), g))

@inline _compose_sf_fc(f, g::FunctionChain) = FunctionChain((g, f))
@inline _compose_sf_fc(f, g::FunctionChain{<:Tuple}) = FunctionChain((g.fs..., f))
@inline _compose_sf_fc(f::F, g::FunctionChain{<:AbstractVector{F}}) where F = FunctionChain(push!(copy(g.fs), f))


_iterate_fs(::Nothing, fs, x) = throw(ArgumentError("Chain of functions must not be an empty iterable"))

function _iterate_fs((f1, itr_state), fs, x)
    y = f1(x)
    next = iterate(fs, itr_state)
    while !isnothing(next)
        f_i, itr_state = next
        y = f_i(y)
        next = iterate(fs, itr_state)
    end
    return y
end

(fc::FunctionChain)(x) = _iterate_fs(iterate(fc.fs), fc.fs, x)

_iterate_fs_withintermediate(::Nothing, fs, x::T) where T = throw(ArgumentError("Chain of functions must not be an empty iterable"))


function _iterate_fs_withintermediate((f1, itr_state), fs, x)
    y = f1(x)
    ys = _similar_empty(fs, typeof(y))
    _sizehint!(ys, Base.IteratorSize(fs), fs)
    ys = _push!!(ys, y)
    next = iterate(fs, itr_state)
    while !isnothing(next)
        f_i, itr_state = next
        y = f_i(y)
        ys = _push!!(ys, y)
        next = iterate(fs, itr_state)
    end
    return ys
end

with_intermediate_results(fc::FunctionChain, x) = _iterate_fs_withintermediate(iterate(fc.fs), fc.fs, x)


function _tuple_fc_exprs(::Type{<:FunctionChain{FS}}) where {FS<:Tuple}
    expr = Expr(:block)
    y_0 = Symbol(:y, 0)
    push!(expr.args, :($y_0 = x))
    idxs = eachindex(FS.parameters)
    for i in idxs
        y_i, x_i = Symbol(:y, i), Symbol(:y, i-1)
        push!(expr.args, :($y_i = fc.fs[$i]($x_i)))
    end
    return expr
end


@inline (fc::FunctionChain{Tuple{}})(x) = x

@inline (fc::FunctionChain{Tuple{F}})(x) where F = fc.fs[1](x)

@generated function (fc::FunctionChain{FS})(x) where {FS<:Tuple}
    expr = _tuple_fc_exprs(fc)
    push!(expr.args, Symbol(:y, last(eachindex(FS.parameters))))
    return expr
end


@inline with_intermediate_results(::FunctionChain{Tuple{}}, x) = ()

@inline with_intermediate_results(fc::FunctionChain{Tuple{F}}, x) where F = (fc.fs[1](x),)

@generated function with_intermediate_results(fc::FunctionChain{FS}, x) where {FS<:Tuple}
    expr = _tuple_fc_exprs(fc)
    push!(expr.args, :(($([Symbol(:y, i) for i in eachindex(FS.parameters)]...),)))
    return expr
end


"""
    fchain()
    fchain(fs)

Construct a function chain of functions `fs`.

Typically returns a [`FunctionChain`](@ref), but may be specialized to return
other function chain types for specific types of functions.

`fs` must be iterable, it may be a tuple, vector, generator, etc.
`fchain(fs)(x)` will apply the functions in `fs` in order of iteration.
"""
function fchain end
export fchain

@inline fchain() = FunctionChain(())

@inline fchain(fs::FS) where FS = _fchain_onearg(fs, Val(static_hasmethod(iterate, Tuple{FS})))
_fchain_onearg(fs::FS, ::Val{true}) where FS = FunctionChain(fs)
_fchain_onearg(f::F, ::Val{false}) where F = FunctionChain(_typed_funcs_tuple(f))

@inline fchain(fs::Tuple{Vararg{Function, N}}) where N = FunctionChain(fs)

function fchain(fs::Tuple{Vararg{Any, N}}) where N
    Base.depwarn("fchain(fs::Tuple) with fs elements not of type Function is deprecated due to possible type instabilities, use `fchain(fs...)` instead.", :fchain)
    fchain(fs...)
end

@inline fchain(fs::Vararg{Any,N}) where N = FunctionChain(_typed_funcs_tuple(fs...))

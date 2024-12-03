# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    with_intermediate_results(f, x)

Apply multi-step function `f` to `x` and return a collection that contains
the intermediate results and the final result.
"""
function with_intermediate_results end
export with_intermediate_results

with_intermediate_results(f, x) = (f(x),)


@inline @generated function _typed_funcs_tuple(fs::Vararg{Any,N}) where N
    expr = Expr(:tuple)
    for i in 1:N
        if fs[i] <: Type
            push!(expr.args, :(AsFunction{$(fs[i])}(fs[$i])))
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

`(fc::FunctionChain)(x)` applies the functions in the chain in order of
iteration over `fc.fs`.

Use [`fchain`](@ref) to construct function chains instead of using the
constructor `FunctionChain(fs)` directly.
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

# Enables `(; fc...)`
Base.merge(a::NamedTuple, fc::FunctionChain{<:NamedTuple{names}}) where {names} = merge(a, fc.fs)

convert(::Type{FunctionChain}, f::ComposedFunction) = ComposedFunction(_flatten_composed(f))

@inline _flatten_composed(f::F) where {F} = (typed_callable(f),)
@inline _flatten_composed(f::F) where {F<:ComposedFunction} = _typed_funcs_tuple(_flatten_composed(f.inner)..., _flatten_composed(f.outer)...)

@inline Base.:(∘)(f::FunctionChain, g::FunctionChain) = _compose_fc_fc(f, g)
@inline Base.:(∘)(f::FunctionChain, g::ComposedFunction) = _compose_fc_fc(f, FunctionChain(_flatten_composed(g)))
@inline Base.:(∘)(f::ComposedFunction, g::FunctionChain) = _compose_fc_fc(FunctionChain(_flatten_composed(f)), g)
@inline Base.:(∘)(f::FunctionChain, g) = _compose_fc_sf(f, typed_callable(g))
@inline Base.:(∘)(f, g::FunctionChain) = _compose_sf_fc(typed_callable(f), g)
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


function _tuple_fc_exprs(::Type{FS}) where {FS<:Tuple}
    expr = Expr(:block)
    y_0 = Symbol(:y, 0)
    push!(expr.args, :($y_0 = x))
    idxs = eachindex(FS.parameters)
    for i in idxs
        y_i, x_i = Symbol(:y, i), Symbol(:y, i-1)
        push!(expr.args, :($y_i = fs[$i]($x_i)))
    end
    return expr
end


@inline (fc::FunctionChain{Tuple{}})(x) = x

@inline (fc::FunctionChain{Tuple{F}})(x) where F = fc.fs[1](x)

@inline (fc::FunctionChain{FS})(x) where {FS<:Tuple} = _tuple_fc_apply(fc.fs, x)

@inline (fc::FunctionChain{<:NamedTuple})(x) = _tuple_fc_apply(values(fc.fs), x)

@generated function _tuple_fc_apply(fs::FS, x) where {FS<:Tuple}
    expr = _tuple_fc_exprs(fs)
    push!(expr.args, Symbol(:y, last(eachindex(FS.parameters))))
    return expr
end


@inline with_intermediate_results(::FunctionChain{Tuple{}}, x) = ()

@inline with_intermediate_results(fc::FunctionChain{Tuple{F}}, x) where F = (fc.fs[1](x),)

@inline with_intermediate_results(fc::FunctionChain{FS}, x) where {FS<:Tuple} = _tuple_fc_with_apply_interm(fc.fs, x)

@inline function with_intermediate_results(fc::FunctionChain{<:NamedTuple{names}}, x) where {names}
    return NamedTuple{names}(_tuple_fc_with_apply_interm(values(fc.fs), x))
end

@generated function _tuple_fc_with_apply_interm(fs::FS, x) where {FS<:Tuple}
    expr = _tuple_fc_exprs(fs)
    push!(expr.args, :(($([Symbol(:y, i) for i in eachindex(FS.parameters)]...),)))
    return expr
end


"""
    fchain()
    fchain(fs)
    fchain(fs...)

Construct a function chain of functions `fs`.

Typically returns a [`FunctionChain`](@ref), but may be specialized to return
other function chain types for specific types of functions.

`fs` must be iterable, it may be a tuple, vector, generator, etc.
`fchain(fs)(x)` will apply the functions in `fs` in order of iteration.

The resulting function chain supports [`with_intermediate_results`](@ref), and
also supports `InverseFunctions.inverse` and/or
`ChangesOfVariables.with_logabsdet_jacobian` if all functions in the chain do so.
"""
function fchain end
export fchain

@inline fchain(;fs...) = _fchain_vararg_impl(values(fs))
@inline _fchain_vararg_impl(::NamedTuple{()}) = fchain(())
@inline _fchain_vararg_impl(fs::NamedTuple) = fchain(fs)

@inline fchain(fs::FS) where FS = _fchain_onearg(fs, Val(static_hasmethod(iterate, Tuple{FS})))
_fchain_onearg(fs::FS, ::Val{true}) where FS = FunctionChain(fs)
_fchain_onearg(f::F, ::Val{false}) where F = FunctionChain(_typed_funcs_tuple(f))

@inline fchain(fs::Tuple{Vararg{Function, N}}) where N = FunctionChain(fs)

function fchain(fs::Tuple{Vararg{Any, N}}) where N
    @info "DEBUG" typeof(fs)
    Base.depwarn("fchain(fs::Tuple) with fs elements not of type Function is deprecated due to possible type instabilities, use `fchain(fs...)` instead.", :fchain)
    fchain(fs...)
end

@inline fchain(fs::Vararg{Any,N}) where N = FunctionChain(_typed_funcs_tuple(fs...))

@inline fchain(fs::NamedTuple{names,<:Tuple{Vararg{Function}}}) where names = FunctionChain(fs)

function fchain(@nospecialize(fs::NamedTuple))
    throw(ArgumentError("Do not use fchain(fs::NamedTuple) or fchain(;fs...) with fs elements not of type Function, due to type instability."))
end


# InverseFunctions support

_reverse(fs) = reverse(fs)
_reverse(fs::Base.Iterators.Repeated) = fs
_reverse(fs::Base.Iterators.Take{<:Base.Iterators.Repeated}) = fs

_is_noinverse(f) = f isa NoInverse
_contains_noinverse(fs) = Val(any(_is_noinverse, fs))
@generated function _contains_noinverse(fs::Tuple)
    result = any(x -> x <: NoInverse, fs.parameters)
    :(Val($result))
end
_contains_noinverse(fs::NamedTuple) = _contains_noinverse(values(fs))
_contains_noinverse(::AbstractArray) = Val(false)
_contains_noinverse(::AbstractArray{<:NoInverse}) = Val(true)
_contains_noinverse(fs::AbstractArray{>:NoInverse}) = Val(any(_is_noinverse, fs))

function InverseFunctions.inverse(fc::FunctionChain)
    inv_fs = map(inverse, _reverse(fc.fs))
    if _contains_noinverse(inv_fs) isa Val{true}
        NoInverse(fc)
    else
        FunctionChain(inv_fs)
    end
end

function InverseFunctions.inverse(fc::FunctionChain{<:Iterators.Take{<:Iterators.Repeated}})
    fs = fc.fs
    return fchain(Iterators.repeated(inverse(fs.xs.x), fs.n))
end

InverseFunctions.inverse(fc::FunctionChain{<:Base.Generator}) = NoInverse(fc)

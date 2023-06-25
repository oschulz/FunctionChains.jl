# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    with_intermediate_results(f, x)

Apply multi-step function `f` to `x` and return a collection that contains
the intermediate results and the final result.
"""
function with_intermediate_results end
export with_intermediate_results

with_intermediate_results(f, x) = (f(x),)


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


function Base.show(io::IO, m::MIME"text/plain", fc::FunctionChain)
    print(io, "FunctionChain(")
    show(io, m, fc.fs)
    print(io, ")")
end

function Base.show(io::IO, fc::FunctionChain)
    print(io, "FunctionChain(")
    show(io, fc.fs)
    print(io, ")")
end


convert(::Type{FunctionChain}, f::ComposedFunction) = ComposedFunction(_flatten_composed(f))

@inline _flatten_composed(f::F) where {F} = (f,)
@inline _flatten_composed(f::F) where {F<:ComposedFunction} = (_flatten_composed(f.inner)..., _flatten_composed(f.outer)...)

@inline Base.:(∘)(f::FunctionChain, g::FunctionChain) = FunctionChain((g, f))
@inline Base.:(∘)(f::FunctionChain, g) = FunctionChain((_flatten_composed(g)..., f))
@inline Base.:(∘)(f, g::FunctionChain) = FunctionChain((g, _flatten_composed(f)...))
@inline Base.:(∘)(f::FunctionChain, ::typeof(identity)) = f
@inline Base.:(∘)(::typeof(identity), g::FunctionChain) = g

@inline Base.:(∘)(f::FunctionChain{<:Tuple}, g::FunctionChain{<:Tuple}) = FunctionChain((g.fs..., f.fs...))
@inline Base.:(∘)(f::FunctionChain{<:Tuple}, g) = FunctionChain((_flatten_composed(g)..., f.fs...))
@inline Base.:(∘)(f, g::FunctionChain{<:Tuple}) = FunctionChain((g.fs..., _flatten_composed(f)...))

@inline Base.:(∘)(f::FunctionChain{<:AbstractVector{F}}, g::FunctionChain{<:AbstractVector{F}}) where F = FunctionChain(vcat(g.fs, f.fs))
@inline Base.:(∘)(f::FunctionChain{<:AbstractVector{F}}, g::F) where F = FunctionChain(pushfirst!(copy(f.fs), _flatten_composed(g)))
#!!! @inline Base.:(∘)(f::F, g::FunctionChain{<:AbstractVector{F}}), where F = FunctionChain(push!(copy(g.fs), _flatten_composed(f)))


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

fchain() = FunctionChain(())
fchain(fs) = FunctionChain(fs)

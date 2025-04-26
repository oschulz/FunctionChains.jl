# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    fchainfs(fc)

Get the component functions of a function chain or composed function `fc`, in
order of function execution.

See [`fchain`](@ref) for details.
"""
function fchainfs end
export fchainfs

@inline fchainfs(f::F) where {F} = (typed_callable(f),)
@inline fchainfs(f::F) where {F<:ComposedFunction} = _typed_funcs_tuple(fchainfs(f.inner)..., fchainfs(f.outer)...)


"""
    with_intermediate_results(f, x)

Apply multi-step function `f` to `x` and return a collection that contains
the intermediate results and the final result.
"""
function with_intermediate_results end
export with_intermediate_results

with_intermediate_results(f, x) = (f(x),)


_typed_funcs_tuple(fs::Vararg{Function,N}) where N = (fs...,)

@inline @generated function _typed_funcs_tuple(fs::Vararg{Any,N}) where N
    expr = Expr(:tuple)
    for i in 1:N
        if fs[i] <: Type
            push!(expr.args, :(AsFunction{$(fs[i])}(fs[$i])))
        # Future option (breaking) - flatten FunctionChains:
        # elseif fs[i] <: FunctionChain{<:Tuple}
        #    push!(expr.args, :(fs[$i]._fs...))
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

`(fc::FunctionChain)(x)` applies the functions in the chain in order of
iteration over `fc._fs`.

Use [`fchain`](@ref) to construct function chains instead of using the
constructor `FunctionChain(fs)` directly.

Use [`fchainfs(fc)`](@ref) to retrieve the components of a
`FunctionChain` in order of function execution.
"""
struct FunctionChain{FS} <: Function
    _fs::FS
end
export FunctionChain


const _BCastedFC{FS} = Base.Broadcast.BroadcastFunction{<:FunctionChain{FS}}


fchainfs(fc::FunctionChain) = getfield(fc, :_fs)


Base.:(==)(a::FunctionChain, b::FunctionChain) = a._fs == b._fs
Base.isapprox(a::FunctionChain, b::FunctionChain; kwargs...) = _isapprox(a._fs, b._fs; kwargs...)

_isapprox(a, b; kwargs...) = isapprox(a, b; kwargs...)
_isapprox(a::Tuple{Vararg{Any,N}}, b::Tuple{Vararg{Any,N}}; kwargs...) where N = all(map((a, b) ->_isapprox(a, b; kwargs...), a, b))


function Base.show(io::IO, m::MIME"text/plain", fc::FunctionChain)
    print(io, "fchain(")
    show(io, m, fc._fs)
    print(io, ")")
end

function Base.show(io::IO, m::MIME"text/plain", fc::FunctionChain{<:Tuple})
    print(io, "fchain")
    show(io, m, fc._fs)
end

Base.show(io::IO, fc::FunctionChain) = show(io, MIME"text/plain"(), fc)

Base.length(fc::FunctionChain) = length(fc._fs)

# Enables `(fc...,)` and `[fc...]`
Base.iterate(fc::FunctionChain) = iterate(fc._fs)
Base.iterate(fc::FunctionChain, state) = iterate(fc._fs, state)

# Enables `(; fc...)`
Base.merge(a::NamedTuple, fc::FunctionChain{<:NamedTuple{names}}) where {names} = merge(a, fc._fs)

convert(::Type{FunctionChain}, f::ComposedFunction) = ComposedFunction(fchainfs(f))

@inline Base.:(∘)(f::FunctionChain, g::FunctionChain) = _compose_fc_fc(f, g)
@inline Base.:(∘)(f::FunctionChain, g::ComposedFunction) = _compose_fc_fc(f, FunctionChain(fchainfs(g)))
@inline Base.:(∘)(f::ComposedFunction, g::FunctionChain) = _compose_fc_fc(FunctionChain(fchainfs(f)), g)
@inline Base.:(∘)(f::FunctionChain, g) = _compose_fc_sf(f, typed_callable(g))
@inline Base.:(∘)(f, g::FunctionChain) = _compose_sf_fc(typed_callable(f), g)
@inline Base.:(∘)(f::FunctionChain, ::typeof(identity)) = f
@inline Base.:(∘)(::typeof(identity), g::FunctionChain) = g

@inline _compose_fc_fc(f::FunctionChain, g::FunctionChain) = FunctionChain((g, f))
@inline _compose_fc_fc(f::FunctionChain{<:Tuple}, g::FunctionChain) = FunctionChain((g, f._fs...))
@inline _compose_fc_fc(f::FunctionChain, g::FunctionChain{<:Tuple}) = FunctionChain((g._fs..., f))
@inline _compose_fc_fc(f::FunctionChain{<:Tuple}, g::FunctionChain{<:Tuple}) = FunctionChain((g._fs..., f._fs...))
@inline _compose_fc_fc(f::FunctionChain{<:AbstractVector{F}}, g::FunctionChain{<:AbstractVector{F}}) where F = FunctionChain(vcat(g._fs, f._fs))

@inline function _compose_fc_fc(
    f::FunctionChain{<:NamedTuple{names_f}}, g::FunctionChain{<:NamedTuple{names_g}}
) where {names_f, names_g}
    _compose_fc_fc_nt(f, g, _no_names_overlap(Val(names_f), Val(names_g)))
end

@inline @generated function _no_names_overlap(::Val{names_a}, ::Val{names_b}) where {names_a, names_b}
    all_names = sort([names_a..., names_b...])
    length(unique(all_names)) == length(all_names) ? :(Val(true)) : :(Val(false))
end

@inline _compose_fc_fc_nt(f, g, ::Val{true}) = FunctionChain((; g._fs..., f._fs...))
@inline _compose_fc_fc_nt(f, g, ::Val{false}) = FunctionChain((g, f))


@inline _compose_fc_sf(f::FunctionChain, g) = FunctionChain((g, f))
@inline _compose_fc_sf(f::FunctionChain{<:Tuple}, g) = FunctionChain((g, f._fs...))
@inline _compose_fc_sf(f::FunctionChain{<:AbstractVector{F}}, g::F) where F = FunctionChain(pushfirst!(copy(f._fs), g))

@inline _compose_sf_fc(f, g::FunctionChain) = FunctionChain((g, f))
@inline _compose_sf_fc(f, g::FunctionChain{<:Tuple}) = FunctionChain((g._fs..., f))
@inline _compose_sf_fc(f::F, g::FunctionChain{<:AbstractVector{F}}) where F = FunctionChain(push!(copy(g._fs), f))


_iterate_fs(::Nothing, fs, x, f_apply) = throw(ArgumentError("Chain of functions must not be an empty iterable"))

function _iterate_fs((f1, itr_state), fs, x, f_apply)
    y = f_apply(f1, x)
    next = iterate(fs, itr_state)
    while !isnothing(next)
        f_i, itr_state = next
        y = f_apply(f_i, y)
        next = iterate(fs, itr_state)
    end
    return y
end

(fc::FunctionChain)(x) = _iterate_fs(iterate(fc._fs), fc._fs, x, applyf)

@inline Base.broadcasted(fc::FunctionChain, x) = _bc_fc(fc, x)

_bc_fc(fc::FunctionChain, x) = _iterate_fs(iterate(fc._fs), fc._fs, x, broadcast)


_iterate_fs_withintermediate(::Nothing, fs, x::T) where T = throw(ArgumentError("Chain of functions must not be an empty iterable"))

function _iterate_fs_withintermediate((f1, itr_state), fs, x, f_apply)
    y = f_apply(f1, x)
    ys = _similar_empty(fs, typeof(y))
    _sizehint!(ys, Base.IteratorSize(fs), fs)
    ys = _push!!(ys, y)
    next = iterate(fs, itr_state)
    while !isnothing(next)
        f_i, itr_state = next
        y = f_apply(f_i, y)
        ys = _push!!(ys, y)
        next = iterate(fs, itr_state)
    end
    return ys
end

function with_intermediate_results(fc::FunctionChain, x)
    fs = fchainfs(fc)
    return _iterate_fs_withintermediate(iterate(fs), fs, x, applyf)
end

function with_intermediate_results(bfc::_BCastedFC, x)
    fs = fchainfs(bfc.f)
    return _iterate_fs_withintermediate(iterate(fs), fs, x, broadcast)
end


function _fc_fs_tpl_expr(
    n::Integer;
    f_apply::Union{Symbol,Expr} = :applyf, postproc::Union{Symbol,Expr} = :identity,
    return_intermediates::Bool = false
)
    if n > 0
        expr = Expr(:block)
        y_syms = Symbol.(:y, 1:n)
        for i in 1:n
            y_i, x_i = y_syms[i], (i > 1 ? y_syms[i-1] : :x)
            if f_apply == :applyf
                push!(expr.args, :($y_i = fs[$i]($x_i)))
            else
                push!(expr.args, :($y_i = ($f_apply)(fs[$i], $x_i)))
            end
        end
        results = if postproc == :identity
            y_syms
        else
            map(ysym -> :($postproc($ysym)), y_syms)
        end
        if return_intermediates
            push!(expr.args, :(return ($(results...),)))
        else
            push!(expr.args, :(return $(last(results))))
        end
        return expr
    else
        if return_intermediates
            return :(return (x,))
        else
            return :(return x)
        end
    end
    
    return expr
end


@inline (fc::FunctionChain{Tuple{}})(x) = x

@inline (fc::FunctionChain{Tuple{F}})(x) where F = fc._fs[1](x)

@inline (fc::FunctionChain{FS})(x) where {FS<:Tuple} = _apply_fc_fs_tpl(fc._fs, x)

@inline (fc::FunctionChain{<:NamedTuple})(x) = _apply_fc_fs_tpl(values(fc._fs), x)

@generated function _apply_fc_fs_tpl(fs::FS, x) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _fc_fs_tpl_expr(n)
end


@inline _bc_fc(fc::FunctionChain{FS}, x) where {FS<:Tuple} = _bc_fc_fs_tpl(fc._fs, x)

@inline _bc_fc(fc::FunctionChain{<:NamedTuple}, x) = _bc_fc_fs_tpl(values(fc._fs), x)

@generated function _bc_fc_fs_tpl(fs::FS, x) where {FS<:Tuple}
    n = length(FS.parameters)
    return _fc_fs_tpl_expr(n, f_apply = :(Base.broadcasted), postproc = :identity, return_intermediates = false)
end


@inline with_intermediate_results(::FunctionChain{Tuple{}}, x) = (x,)

@inline with_intermediate_results(fc::FunctionChain{Tuple{F}}, x) where F = (fc._fs[1](x),)

@inline with_intermediate_results(fc::FunctionChain{FS}, x) where {FS<:Tuple} = _apply_interm_fc_fs_tpl(fc._fs, x)

@inline function with_intermediate_results(fc::FunctionChain{<:NamedTuple{names}}, x) where {names}
    return NamedTuple{names}(_apply_interm_fc_fs_tpl(values(fc._fs), x))
end

@generated function _apply_interm_fc_fs_tpl(fs::FS, x) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _fc_fs_tpl_expr(n, return_intermediates = true)
end


with_intermediate_results(bfc::_BCastedFC{FS}, x) where {FS<:Tuple} = _bc_interm_fc_fs_tpl(bfc.f._fs, x)

@inline function with_intermediate_results(bfc::_BCastedFC{<:NamedTuple{names}}, x) where {names}
    return NamedTuple{names}(_bc_interm_fc_fs_tpl(values(bfc.f._fs), x))
end

@generated function _bc_interm_fc_fs_tpl(fs::FS, x) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _fc_fs_tpl_expr(n, f_apply = :(Base.broadcast), postproc = :identity, return_intermediates = true)
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

Use [`fchainfs(fc)`](@ref) to retrieve `fs`.
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
_reverse(fs::Repeated) = fs
_reverse(fs::Take{<:Repeated}) = fs

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
    inv_fs = map(inverse, _reverse(fc._fs))
    if _contains_noinverse(inv_fs) isa Val{true}
        NoInverse(fc)
    else
        FunctionChain(inv_fs)
    end
end

function InverseFunctions.inverse(fc::FunctionChain{<:Iterators.Take{<:Iterators.Repeated}})
    fs = fc._fs
    return fchain(Iterators.repeated(inverse(fs.xs.x), fs.n))
end

InverseFunctions.inverse(fc::FunctionChain{<:Base.Generator}) = NoInverse(fc)



"""
    ffchain(f, g, hs...)
    ffchain() = identity
    ffchain(f) = f
    ffchain(::Type{F}) where F = FunctionChains.AsFunction{Type{F}}(F)

Similar to [`fchain((f, g, hs...))`](@ref), but flattens arguments of type
`ComposedFunction` and merges [`FunctionChain`](@ref) arguments.

Tries to remove superfluous `identity` functions and to return a simple
function instead of a `FunctionChain` if possible.

Behaves like `ffcomp(hs..., g, f)` (see [`ffcomp`](@ref)).
"""
function ffchain end
export ffchain

@inline ffchain() = identity
@inline ffchain(f) = f
@inline ffchain(::Type{F}) where F = FunctionChains.AsFunction{Type{F}}(F)
@inline ffchain(f::ComposedFunction) = _flat_fs_postproc(_flat_fs(f))

@inline _flat_fs(f::F) where F = (f,)
@inline _flat_fs(::typeof(identity)) = ()
@inline _flat_fs(::Type{F}) where F = (FunctionChains.AsFunction{Type{F}}(F),)
@inline _flat_fs(f::ComposedFunction) = (_flat_fs(f.inner)..., _flat_fs(f.outer)...)
@inline _flat_fs(f::FunctionChain{<:Tuple}) = fchainfs(f)

@inline _flat_fs_postproc(::Tuple{}) = identity
@inline _flat_fs_postproc(fs::Tuple{F}) where F = fs[1]
@inline _flat_fs_postproc(fs::Tuple) = FunctionChain(fs)

@inline @generated function ffchain(fs::Vararg{Any,N}) where N
    expr = Expr(:tuple)
    for i in 1:N
        if !(fs[i] <: typeof(identity))
            push!(expr.args, :(_flat_fs(fs[$i])...))
        end
    end
    if length(expr.args) == 0
        return :(identity)
    elseif length(expr.args) == 1
        only_fs = only(only(expr.args).args)
        return :(_flat_fs_postproc($only_fs))
    else
        return :(FunctionChain($(expr)))
    end
end

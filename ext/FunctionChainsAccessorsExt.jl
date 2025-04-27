# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsAccessorsExt

import Accessors
using Accessors: set

using InverseFunctions: inverse, NoInverse

using Base: IteratorSize

using FunctionChains
using FunctionChains: _fc_fs_tpl_expr, _BCastedFC, _iterate_fs_withintermediate
using FunctionChains: _check_fp_sizes


# FunctionChain ==============================================================

Accessors.set(x, fc::FunctionChain, z) = _set_maybe_withinverse(fc, inverse(fc), x, z)
Accessors.set(x, fc::_BCastedFC, z) = _set_maybe_withinverse(fc, inverse(fc), x, z)

_set_maybe_withinverse(@nospecialize(fc::FunctionChain), inv_fc, @nospecialize(x), z) = inv_fc(z)
_set_maybe_withinverse(@nospecialize(fc::_BCastedFC), inv_fc, @nospecialize(x), z) = inv_fc(z)


@inline _set_maybe_withinverse(fc::FunctionChain{FS}, ::NoInverse, x, z) where {FS<:Tuple} = _set_fc_fs_tpl(fc._fs, x, z)
@inline _set_maybe_withinverse(fc::FunctionChain{<:NamedTuple}, ::NoInverse, x, z) = _set_fc_fs_tpl(values(fc._fs), x, z)

@generated function _set_fc_fs_tpl(fs::FS, x, z) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _fc_fs_tpl_expr(n, f_apply = :applyf, postproc=:identity, reduction=:accessors_set)
end

@inline _set_maybe_withinverse(bfc::_BCastedFC{FS}, ::NoInverse, x, z) where {FS<:Tuple} = _set_bc_fc_fs_tpl(bfc.f._fs, x, z)
@inline _set_maybe_withinverse(bfc::_BCastedFC{<:NamedTuple}, ::NoInverse, x, z) = _set_bc_fc_fs_tpl(values(bfc.f._fs), x, z)

@generated function _set_bc_fc_fs_tpl(fs::FS, x, z) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _fc_fs_tpl_expr(n, f_apply = :(Base.broadcast), postproc=:identity, reduction=:accessors_set)
end



@inline _set_maybe_withinverse(fc::FunctionChain, ::NoInverse, x, z) = _set_fc_fs_iterable(fc._fs, x, z, applyf, identity)
@inline _set_maybe_withinverse(fc::_BCastedFC, ::NoInverse, x, z) = _set_fc_fs_iterable(fc.f._fs, x, z, broadcast, fbcast)

function _set_fc_fs_iterable(fs, x, z, f_apply, f_wrap)
    n = length(fs)
    if n < 1
        throw(ArgumentError("Chain of functions must not be an empty iterable"))
    elseif n == 1
        return set(x, first(fs), z)
    else
        short_fs = _without_last(fs)
        ys = _iterate_fs_withintermediate(iterate(short_fs), short_fs, x, f_apply)

        rev_fs = Iterators.reverse(fs)
        rev_ys = Iterators.reverse(ys)

        f_i, istate_f = iterate(rev_fs)
        y_i, istate_y = iterate(rev_ys)

        wrapped_f_i = f_wrap(f_i)
        result = set(y_i, wrapped_f_i, z)

        next_f = iterate(rev_fs, istate_f)
        next_y = iterate(rev_ys, istate_y)
        while !isnothing(next_f)
            f_i, istate_f = next_f
            wrapped_f_i = f_wrap(f_i)

            if !isnothing(next_y)
                y_i, istate_y = next_y
                result = set(y_i, wrapped_f_i, result)
                next_y = iterate(rev_ys, istate_y)
            else
                result = set(x, wrapped_f_i, result)
            end
            next_f = iterate(rev_fs, istate_f)
        end
        return result
    end
end

_without_last(fs::AbstractArray) = view(fs, firstindex(fs):lastindex(fs)-1)
_without_last(fs) = Iterators.take(fs, length(fs)-1)


# FCartProd ==============================================================

Accessors.set(x, fp::FCartProd, z) = _fp_set(fp._fs, x, z)


function _fp_set(fs::Tuple{Vararg{Any,N}}, x::Tuple{Vararg{Any,N}}, z::Tuple{Vararg{Any,N}}) where N
    map(set, x, fs, z)
end

_fp_set(::Tuple{Vararg{typeof(identity),N}}, @nospecialize(x::Tuple{Vararg{Any,N}}), z::Tuple{Vararg{Any,N}}) where N = z

function _fp_set(@nospecialize(fs::Tuple), @nospecialize(x::Tuple), @nospecialize(z::Tuple))
    throw(ArgumentError("Can't Accessors.set a Tuple of length $(length(x)) via FCartProd over Tuple of length $(length(fs)) with a Tuple of length $(length(z))."))
end


function _fp_set(fs::NamedTuple{names}, x::NamedTuple{names}, z::NamedTuple{names}) where names
    NamedTuple{names}(map(set, values(x), values(fs), values(z)))
end

_fp_set(::NamedTuple{names,<:Tuple{Vararg{typeof(identity)}}}, x::NamedTuple{names}, @nospecialize(z::NamedTuple{names})) where names = z

function _fp_set(@nospecialize(fs::NamedTuple), @nospecialize(x::NamedTuple), @nospecialize(z::NamedTuple))
    throw(ArgumentError("Can't Accessors.set a NamedTuple with names $(propertynames(x)) via FCartProd over NamedTuple with names $(propertynames(fs)) with a NamedTuple with names $(propertynames(z))."))
end


Base.@propagate_inbounds function _fp_set(fs, x, z)
    @boundscheck _check_fp_sizes(IteratorSize(fs), IteratorSize(x), fs, x)
    @inbounds result = _fp_set_generic_impl(fs, x, z)
    return result
end

Base.@propagate_inbounds _fp_set_generic_impl(fs, x, z) = set.(x, fs, z)

_fp_set_generic_impl(::AbstractArray{typeof(identity)}, @nospecialize(x), z) = z


end # module FunctionChainsAccessorsExt

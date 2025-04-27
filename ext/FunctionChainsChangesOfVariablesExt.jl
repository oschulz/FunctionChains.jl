# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsChangesOfVariablesExt

using ChangesOfVariables
using ChangesOfVariables: with_logabsdet_jacobian
using FunctionChains

using FunctionChains: AsFunction
using FunctionChains: _BCastedFC, _check_fp_sizes

using Base: IteratorSize, HasShape, HasLength

using Core: Typeof


# AsFunction ==============================================================

ChangesOfVariables.with_logabsdet_jacobian(ff::AsFunction, x) = with_logabsdet_jacobian(ff.f, x)


# FunctionChain ==============================================================

@inline ChangesOfVariables.with_logabsdet_jacobian(fc::FunctionChain, x) = _withladj_fc(fc, x)

_iterate_fs_withladj(::Nothing, fs, x, f_wrap) = throw(ArgumentError("Chain of functions must not be an empty iterable"))

function _iterate_fs_withladj((f1, itr_state), fs, x, f_wrap)
    y_ladj = with_logabsdet_jacobian(f_wrap(f1), x)
    y_ladj isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{FunctionChain{Typeof(fs)},Typeof(x)}()
    y, ladj = y_ladj
    next = iterate(fs, itr_state)
    while !isnothing(next)
        f_i, itr_state = next
        y_ladj_i = with_logabsdet_jacobian(f_wrap(f_i), y)
        y_ladj_i isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{FunctionChain{Typeof(fs)},Typeof(x)}()
        y, ladj_i = y_ladj_i
        ladj += ladj_i
        next = iterate(fs, itr_state)
    end
    return y, ladj
end


function _withladj_fc(fc::FunctionChain, x)
    fs = fchainfs(fc)
    return _iterate_fs_withladj(iterate(fs), fs, x, identity)
end

function _withladj_fc(bfc::_BCastedFC, x)
    fs = fchainfs(bfc.f)
    return _iterate_fs_withladj(iterate(fs), fs, x, fbcast)
end


@inline _withladj_fc(::FunctionChain{Tuple{}}, x) = with_logabsdet_jacobian(identity, x)

@inline _withladj_fc(fc::FunctionChain{Tuple{F}}, x) where F = with_logabsdet_jacobian(fc._fs[1], x)

@inline _withladj_fc(fc::FunctionChain{FS}, x) where {FS<:Tuple} = _withladj_fc_fs_tpl(fc._fs, x, NoLogAbsDetJacobian{Typeof(fc),Typeof(x)}())

@inline _withladj_fc(fc::FunctionChain{<:NamedTuple}, x) = _withladj_fc_fs_tpl(values(fc._fs), x, NoLogAbsDetJacobian{Typeof(fc),Typeof(x)}())

@generated function _withladj_fc_fs_tpl(fs::FS, x, no_ladj) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _withladj_fc_fs_tpl_expr(n, f_wrap = :identity)
end


@inline _withladj_fc(bfc::_BCastedFC{FS}, x) where {FS<:Tuple} = _withladj_fc_fs_tpl(bfc.f._fs, x, NoLogAbsDetJacobian{Typeof(fc),Typeof(x)}())

@inline function _withladj_fc(bfc::_BCastedFC{<:NamedTuple{names}}, x) where {names}
    return NamedTuple{names}(_withladj_fc_fs_tpl(values(bfc.f._fs), x, NoLogAbsDetJacobian{Typeof(fc),Typeof(x)}()))
end

@generated function _withladj_bcfc_fs_tpl(fs::FS, x, no_ladj) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _withladj_fc_fs_tpl_expr(n, f_wrap = :fbcast)
end


function _withladj_fc_fs_tpl_expr(n::Integer; f_wrap::Union{Symbol,Expr} = :identity)
    if n == 0
        return :(with_logabsdet_jacobian(identity, x))
    else
        expr = Expr(:block)
        y_syms = Symbol.(:y, 1:n)
        ladj_syms = Symbol.(:ladj, 1:n)
        y_ladj_syms = Symbol.(:y_ladj, 1:n)
        for i = 1:n
            sym_input = i > 1 ? y_syms[i-1] : :x
            f_expr = (f_wrap == :identity) ? :(fs[$i]) : :($(f_wrap)(fs[$i]))
            push!(expr.args, :($(y_ladj_syms[i]) = with_logabsdet_jacobian($f_expr, $sym_input)))
            push!(expr.args, :($(y_ladj_syms[i]) isa NoLogAbsDetJacobian && return no_ladj))
            push!(expr.args, :($(y_syms[i]) = $(y_ladj_syms[i])[1]))
            if i > 1
                push!(expr.args, :($(ladj_syms[i]) = $(y_ladj_syms[i])[2] + $(ladj_syms[i-1])))
            else
                push!(expr.args, :($(ladj_syms[i]) = $(y_ladj_syms[i])[2]))
            end
        end
        push!(expr.args, :(return ($(last(y_syms)), $(last(ladj_syms)))))
        return expr
    end
end



# FCartProd ==============================================================


_is_noladj(r) = r isa NoLogAbsDetJacobian
_contains_noladj(r::NoLogAbsDetJacobian) = Val(true)
_contains_noladj(r) = Val(any(_is_noladj, r))
@generated function _contains_noladj(r::Tuple)
    result = any(x -> x <: NoLogAbsDetJacobian, r.parameters)
    :(Val($result))
end
_contains_noladj(r::NamedTuple) = _contains_noladj(values(r))
_contains_noladj(::AbstractArray) = Val(false)
_contains_noladj(::AbstractArray{<:NoLogAbsDetJacobian}) = Val(true)
_contains_noladj(r::AbstractArray{>:NoLogAbsDetJacobian}) = Val(any(_is_noladj, r))


# Partially from ChangesOfVariables.jl, modified/improved:

@inline _sum_ladjs(@Base.nospecialize(F::Type), @Base.nospecialize(T::Type), ys_with_ladjs::Tuple{Any,Real}) = ys_with_ladjs
@inline _sum_ladjs(::Type{F}, ::Type{T}, ::NoLogAbsDetJacobian) where {F,T} = NoLogAbsDetJacobian{F,T}()

_get_all_first(x) = map(first, x)
# Use x -> x[2] instead of last, using last causes horrible performance in Zygote here:
_sum_over_second(x) = sum(x -> x[2], x)

function _sum_ladjs(::Type{F}, ::Type{T}, ys_with_ladjs) where {F,T}
    if _contains_noladj(ys_with_ladjs) isa Val{true}
        return NoLogAbsDetJacobian{F,T}()
    else
        return _sum_ladjs_impl(ys_with_ladjs)
    end
end

function _sum_ladjs_impl(ys_with_ladjs)
    y = _get_all_first(ys_with_ladjs)
    ladj = _sum_over_second(ys_with_ladjs)
    return (y, ladj)
end


# Implementation of with_logabsdet_jacobian on FCartProd:

function ChangesOfVariables.with_logabsdet_jacobian(fp::FCartProd{<:Tuple{Vararg{Any,N}}}, x::Tuple{Vararg{Any,N}}) where N
    return _sum_ladjs(typeof(fp), typeof(x), map(with_logabsdet_jacobian, fp._fs, x))
end

function ChangesOfVariables.with_logabsdet_jacobian(::FCartProd{<:Tuple{Vararg{typeof(identity),N}}}, x::Tuple{Vararg{Any,N}}) where N
    return with_logabsdet_jacobian(identity, x)
end

function ChangesOfVariables.with_logabsdet_jacobian(@nospecialize(fs::FCartProd{<:Tuple}), @nospecialize(x::Tuple))
    throw(ArgumentError("Can't apply FCartProd over Tuple of length $(length(fp._fs)) to Tuple of length $(length(x))."))
end


function ChangesOfVariables.with_logabsdet_jacobian(fp::FCartProd{<:NamedTuple{names}}, x::NamedTuple{names}) where names
    return _sum_ladjs(typeof(fp), typeof(x), map(with_logabsdet_jacobian, fp._fs, x))
end

function ChangesOfVariables.with_logabsdet_jacobian(::FCartProd{<:NamedTuple{names,<:Tuple{Vararg{typeof(identity)}}}}, x::NamedTuple{names}) where names
    return with_logabsdet_jacobian(identity, x)
end

function ChangesOfVariables.with_logabsdet_jacobian(@nospecialize(fs::FCartProd{<:NamedTuple}), @nospecialize(x::NamedTuple))
    throw(ArgumentError("Can't apply FCartProd over NamedTuple with names $(propertynames(fp._fs)) to NamedTuple with names $(propertynames(x))."))
end


Base.@propagate_inbounds function ChangesOfVariables.with_logabsdet_jacobian(fp::FCartProd, x)
    fs = fp._fs
    @boundscheck _check_fp_sizes(IteratorSize(fs), IteratorSize(x), fs, x)
    @inbounds result = _fp_apply_withladj(fp, x)
    return result
end

Base.@propagate_inbounds function _fp_apply_withladj(fp::FCartProd, x)
    y_with_ladj = with_logabsdet_jacobian.(fp._fs, x)
    _sum_ladjs(typeof(fp), typeof(x), y_with_ladj)
end

_fp_apply_withladj(::FCartProd{<:AbstractArray{typeof(identity)}}, x) = with_logabsdet_jacobian(identity, x)


end # module FunctionChainsChangesOfVariablesExt

# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsChangesOfVariablesExt

using ChangesOfVariables
using FunctionChains

using FunctionChains: _BCastedFC

@inline ChangesOfVariables.with_logabsdet_jacobian(fc::FunctionChain, x) = _withladj_fc(fc, x)


_iterate_fs_withladj(::Nothing, fs, x, f_wrap) = throw(ArgumentError("Chain of functions must not be an empty iterable"))

function _iterate_fs_withladj((f1, itr_state), fs, x, f_wrap)
    y_ladj = with_logabsdet_jacobian(f_wrap(f1), x)
    y_ladj isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{FunctionChain{typeof(fs)},typeof(x)}()
    y, ladj = y_ladj
    next = iterate(fs, itr_state)
    while !isnothing(next)
        f_i, itr_state = next
        y_ladj_i = with_logabsdet_jacobian(f_wrap(f_i), y)
        y_ladj_i isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{FunctionChain{typeof(fs)},typeof(x)}()
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

@inline _withladj_fc(fc::FunctionChain{FS}, x) where {FS<:Tuple} = _withladj_fc_fs_tpl(fc._fs, x, NoLogAbsDetJacobian{typeof(fc),typeof(x)}())

@inline _withladj_fc(fc::FunctionChain{<:NamedTuple}, x) = _withladj_fc_fs_tpl(values(fc._fs), x, NoLogAbsDetJacobian{typeof(fc),typeof(x)}())

@generated function _withladj_fc_fs_tpl(fs::FS, x, no_ladj) where {FS<:Tuple}
    n = length(eachindex(FS.parameters))
    return _withladj_fc_fs_tpl_expr(n, f_wrap = :identity)
end


@inline _withladj_fc(bfc::_BCastedFC{FS}, x) where {FS<:Tuple} = _withladj_fc_fs_tpl(bfc.f._fs, x, NoLogAbsDetJacobian{typeof(fc),typeof(x)}())

@inline function _withladj_fc(bfc::_BCastedFC{<:NamedTuple{names}}, x) where {names}
    return NamedTuple{names}(_withladj_fc_fs_tpl(values(bfc.f._fs), x, NoLogAbsDetJacobian{typeof(fc),typeof(x)}()))
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


end # module FunctionChainsChangesOfVariablesExt

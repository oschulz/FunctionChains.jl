# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsChangesOfVariablesExt

using ChangesOfVariables
using FunctionChains


_iterate_fs_withladj(::Nothing, fs, x) = throw(ArgumentError("Chain of functions must not be an empty iterable"))

function _iterate_fs_withladj((f1, itr_state), fs, x)
    y_ladj = with_logabsdet_jacobian(f1, x)
    y_ladj isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{FunctionChain{typeof(fs)},typeof(x)}()
    y, ladj = y_ladj
    next = iterate(fs, itr_state)
    while !isnothing(next)
        f_i, itr_state = next
        y_ladj_i = with_logabsdet_jacobian(f_i, y)
        y_ladj_i isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{FunctionChain{typeof(fs)},typeof(x)}()
        y, ladj_i = y_ladj_i
        ladj += ladj_i
        next = iterate(fs, itr_state)
    end
    return y, ladj
end


ChangesOfVariables.with_logabsdet_jacobian(fc::FunctionChain, x) = _iterate_fs_withladj(iterate(fc.fs), fc.fs, x)


@inline ChangesOfVariables.with_logabsdet_jacobian(::FunctionChain{Tuple{}}, x) = with_logabsdet_jacobian(identity, x)

@inline ChangesOfVariables.with_logabsdet_jacobian(fc::FunctionChain{Tuple{F}}, x) where F = with_logabsdet_jacobian(fc.fs[1], x)

@generated function ChangesOfVariables.with_logabsdet_jacobian(fc::FunctionChain{FS}, x) where {FS<:Tuple}
    expr = Expr(:block)
    sym_y, sym_ladj, sym_y_ladj, sym_tmp_ladj = gensym(:y), gensym(:ladj), gensym(:y_ladj), gensym(:tmp_ladj)
    push!(expr.args, :($sym_y_ladj = with_logabsdet_jacobian(fc.fs[1], x)))
    push!(expr.args, :($sym_y_ladj isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{typeof(fc),typeof(x)}()))
    push!(expr.args, :(($sym_y, $sym_ladj) = $sym_y_ladj))
    sym_last_y, sym_last_ladj = sym_y, sym_ladj
    for i = 2:length(FS.parameters)
        sym_y, sym_ladj, sym_y_ladj, sym_tmp_ladj = gensym(:y), gensym(:ladj), gensym(:y_ladj), gensym(:tmp_ladj)
        push!(expr.args, :($sym_y_ladj = with_logabsdet_jacobian(fc.fs[$i], $sym_last_y)))
        push!(expr.args, :($sym_y_ladj isa NoLogAbsDetJacobian && return NoLogAbsDetJacobian{typeof(fc),typeof(x)}()))
        push!(expr.args, :(($sym_y, $sym_tmp_ladj) = $sym_y_ladj))
        push!(expr.args, :($sym_ladj = $sym_tmp_ladj + $sym_last_ladj))
        sym_last_y, sym_last_ladj = sym_y, sym_ladj
    end
    push!(expr.args, :(return ($sym_y, $sym_ladj)))

    return expr
end


end # module FunctionChainsChangesOfVariablesExt

# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    fcomp()
    fcomp(fs)
    fcomp(fs...)

Construct a composition of the functions `fs` that is semantically equivalent
to `fs[1] ∘ fs[2] ∘ ...` and `fchain(reverse(fs))`.

The resulting function composition supports
[`with_intermediate_results`](@ref), but note that the intermediate results
are in order of function evaluation, not in the order of the functions in
`fs`.

The composition also supports `InverseFunctions.inverse` and/or
`ChangesOfVariables.with_logabsdet_jacobian` if all functions in the
composition do so.
"""
function fcomp end
export fcomp

@inline fcomp() = fchain(())

@inline fcomp(fs) = fchain(reverse(fs))
@inline fcomp(fs::Tuple{Vararg{Function}}) = fchain(reverse(fs))

@inline fcomp(fs::Vararg{Any}) = fchain(reverse(fs)...)

function fcomp(@nospecialize(fs::Tuple))
    throw(ArgumentError("Do not use fcomp(fs::Tuple) with fs elements not of type Function, due to possible type instabilities, use `fcomp(fs...)` instead."))
end


"""
    ffcomp(f, g, hs...)
    ffcomp() = identity
    ffcomp(f) = f
    ffcomp(::Type{F}) where F = FunctionChains.AsFunction{Type{F}}(F)

Similar to [`fcomp((f, g, hs...))`](@ref), but flattens arguments of type
`ComposedFunction` and merges [`FunctionChain`](@ref) arguments.

Tries to remove superfluous `identity` functions and to return a simple
function instead of a `FunctionChain` if possible.

Behaves like `ffchain(hs..., g, f)` (see [`ffchain`](@ref)).
"""
function ffcomp end
export ffcomp

@inline ffcomp() = identity
@inline ffcomp(f) = f
@inline ffcomp(::Type{F}) where F = FunctionChains.AsFunction{Type{F}}(F)
@inline ffcomp(f::ComposedFunction) = _flat_fs_postproc(_flat_fs(f))

@inline ffcomp(fs::Vararg{Any,N}) where N = ffchain(reverse(fs)...)

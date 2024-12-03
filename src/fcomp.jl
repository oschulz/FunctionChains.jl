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

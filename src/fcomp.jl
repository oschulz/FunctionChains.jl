# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    fcomp()
    fcomp(fs)
    fcomp(fs...)

Construct a composition of the functions `fs` that is functionally equivalent
to `fs[1] ∘ fs[2] ∘ ...` and `fchain(reverse(fs))`.
"""
function fcomp end
export fcomp

@inline fcomp() = fchain(())

@inline fcomp(fs) = fchain(reverse(fs))
@inline fcomp(fs::Tuple{Vararg{Function}}) = fchain(reverse(fs))

@inline fcomp(fs::Vararg{Any}) = fchain(reverse(fs)...)

function fcomp(@nospecialize(fs::Tuple{Vararg{Any}}))
    throw(ArgumentError("Do not use fcomp(fs::Tuple) with fs elements not of type Function, due to possible type instabilities, use `fcomp(fs...)` instead."))
end

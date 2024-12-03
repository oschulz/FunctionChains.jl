# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    frepeat(f, n::Integer)

Create a chain of function `f` repeated `n` times.

See also [`∘̂(f, n)`](@ref) for a version that supports `n <= 0`.

The resulting function chain supports `InverseFunctions.inverse` and/or
`ChangesOfVariables.with_logabsdet_jacobian` if `f` does so.
"""
frepeat(f, n::Integer) = fchain(Iterators.repeated(f, n))
export frepeat


"""
    ∘̂(f, n::Integer)
    f∘̂ n

Create a chain of function `f` repeated `n` times, with support for `n <= 0`.

If `n > 0` then `f∘̂ n` behaves like [`frepeat(f, n)`](@ref). If `n == 0` then
it returns the identity function and if `n < 0` then it behaves like
`frepeat(InverseFunctions.inverse(f), n)`.
"""
@inline function ∘̂(f, n::Integer)
    if n == 0
        return identity
    elseif n > 0
        return frepeat(f, n)
    else
        return frepeat(inverse(f), -n)
    end
end
export ∘̂

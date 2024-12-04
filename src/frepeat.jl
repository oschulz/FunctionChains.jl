# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    frepeat(f, n::Integer)

Create a chain of function `f` repeated `n` times.

See also [`∘̂(f, n)`](@ref) for a version that supports `n <= 0`.

The resulting function chain supports `InverseFunctions.inverse` and/or
`ChangesOfVariables.with_logabsdet_jacobian` if `f` does so.
"""
function frepeat end
export frepeat

function frepeat(f, n::Integer)
    n >= 0 || throw(ArgumentError("frepeat(f, n) requires n to be non-negative"))
    fchain(Iterators.repeated(f, n))
end

function frepeat(f::FunctionChain{<:Take{<:Repeated}}, n::Integer)
    fs = fchainfs(f)
    n_current = length(fs)
    f_orig = fs.xs.x
    return frepeat(f_orig, n_current * n)
end


"""
    ∘̂(f, n::Integer)
    f∘̂ n

Create a chain of function `f` repeated `n` times, with support for `n <= 0`.

The type of the returned function depends on `n`:

* `n == 0`: return `identity`
* `n == 1`: return `f`
* `n > 1`: return `frepeat(f, n)`
* `n == -1`: return `InverseFunctions.inverse(f)`
* `n > -1`: return `frepeat(InverseFunctions.inverse(f), -n)`
"""
@inline function ∘̂(f, n::Integer)
    if n == 0
        return identity
    elseif n == 1
        return f
    elseif n == -1
        return inverse(f)
    elseif n > 0
        return frepeat(f, n)
    else
        return frepeat(inverse(f), -n)
    end
end
export ∘̂

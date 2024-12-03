# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


"""
    fbcast(f)

Return a broadcasted version of the function `f`, so that
`fbcast(f)(A, ...)` is semantically equivalent to `f.(A, ...)`.

`fbcast(f)(A)` is also semantically equivalent to
`fcprod( Fill(f, size(A)) )(A)`.

Typically returns a `Broadcast.BroadcastFunction`, but may be specialized to
return broadcasted implementations depending on the type of `f`. For example,
`fbcast(identity) === identity`.

The resulting broadcasted function supports `InverseFunctions.inverse` and/or
`ChangesOfVariables.with_logabsdet_jacobian` if `f` does so.
"""
function fbcast end
export fbcast

@inline fbcast(f::F) where F = Broadcast.BroadcastFunction{F}(f)
@inline fbcast(::Type{F}) where F = Broadcast.BroadcastFunction{Type{F}}(F)

@inline fbcast(::typeof(identity)) = identity

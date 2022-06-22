# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

import LinearAlgebra, Random
import InverseFunctions, ChangesOfVariables


struct AffineStep{F<:Union{typeof(+),typeof(-),typeof(*),typeof(\)},T}
    a::T
end

AffineStep(f::F, a::T) where {F,T} = AffineStep{F,T}(a)

(f::AffineStep{typeof(+)})(x) = x .+ f.a
(f::AffineStep{typeof(-)})(x) = x .- f.a
(f::AffineStep{typeof(*)})(x) = f.a * x
(f::AffineStep{typeof(\)})(x) = f.a \ x

InverseFunctions.inverse(f::AffineStep{typeof(+)}) = AffineStep(-, f.a)
InverseFunctions.inverse(f::AffineStep{typeof(-)}) = AffineStep(+, f.a)
InverseFunctions.inverse(f::AffineStep{typeof(*)}) = AffineStep(\, f.a)
InverseFunctions.inverse(f::AffineStep{typeof(\)}) = AffineStep(*, f.a)

# Julia v1.8 supports logabsdet(::Number), but older versions don't:
_logabsdet(x::Number) = log(abs(x))
_logabsdet(x::AbstractMatrix) = first(LinearAlgebra.logabsdet(x))

_type_ndof(::Type{<:Real}) = 1
_type_ndof(::Type{<:Complex}) = 2

_mul_ladj(a, x) = _logabsdet(a) * length(eachindex(x)) / length(axes(a,1)) * _type_ndof(eltype(x))

_realtype(::Type{T}) where {T<:Real} = T
_realtype(::Type{Complex{T}}) where {T<:Real} = T

const RCNumber = Union{Real,Complex}
#const RCVector = AbstractArray{<:RCNumber,1}
#const RCMatrix = AbstractArray{<:RCNumber,2}
#const RCArray = AbstractArray{<:RCNumber}

ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(+)}, x) = f(x), zero(_realtype(eltype(x)))
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(-)}, x) = f(x), zero(_realtype(eltype(x)))

ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:Real}, x::Union{RCNumber,Array{<:RCNumber}}) = f(x), _mul_ladj(f.a, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:Complex}, x::Union{Complex,Array{<:Complex}}) = f(x), _mul_ladj(f.a, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:AbstractMatrix{<:Real}}, x::Union{AbstractVector{<:RCNumber},AbstractMatrix{<:RCNumber}}) = f(x), _mul_ladj(f.a, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(*),<:AbstractMatrix{<:Complex}}, x::Union{AbstractVector{<:Complex},AbstractMatrix{<:Complex}}) = f(x), _mul_ladj(f.a, x)

ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:Real}, x::Union{RCNumber,Array{<:RCNumber}}) = f(x), - _mul_ladj(f.a, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:Complex}, x::Union{Complex,Array{<:Complex}}) = f(x), - _mul_ladj(f.a, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:AbstractMatrix{<:Real}}, x::Union{AbstractVector{<:RCNumber},AbstractMatrix{<:RCNumber}}) = f(x), - _mul_ladj(f.a, x)
ChangesOfVariables.with_logabsdet_jacobian(f::AffineStep{typeof(\),<:AbstractMatrix{<:Complex}}, x::Union{AbstractVector{<:Complex},AbstractMatrix{<:Complex}}) = f(x), - _mul_ladj(f.a, x)



struct LeapfrogIntegrator{AF,T}
    A::AF
    Δt::T
end

function (f::LeapfrogIntegrator)(state)
    A, Δt = f.A, f.Δt
    x, v = state.x, state.v  # x[i], v_[i-1/2]
    v_new = v + A(x) * Δt
    x_new = x + v_new*Δt
    return (x = x_new, v = v_new)
end

struct InvLeapfrogIntegrator{AF,T}
    A::AF
    Δt::T
end

function (f::InvLeapfrogIntegrator)(state)
    A, Δt = f.A, f.Δt
    x_new, v_new = state.x, state.v  # x[i+1], v_[i+1/2]
    x = x_new - v_new*Δt
    v = v_new - A(x) * Δt
    return (x = x, v = v)
end

InverseFunctions.inverse(f::LeapfrogIntegrator) = InvLeapfrogIntegrator(f.A, f.Δt)
InverseFunctions.inverse(f::InvLeapfrogIntegrator) = LeapfrogIntegrator(f.A, f.Δt)

ChangesOfVariables.with_logabsdet_jacobian(f::LeapfrogIntegrator, state) = f(state), zero(promote_type(eltype(state.x), eltype(state.v)))
ChangesOfVariables.with_logabsdet_jacobian(f::InvLeapfrogIntegrator, state) = f(state), zero(promote_type(eltype(state.x), eltype(state.v)))



struct MCStep{F,RNG<:Random.AbstractRNG}
    f::F
    rng::RNG
end

MCStep(f) = MCStep(f, Random.default_rng())

(mcstep::MCStep)(x) = mcstep.f(mcstep.rng, x)

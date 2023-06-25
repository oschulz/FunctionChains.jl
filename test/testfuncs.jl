# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

import LinearAlgebra, Random
import InverseFunctions, ChangesOfVariables


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

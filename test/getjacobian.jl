# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

if !isdefined(Main, :getjacobian)

import ForwardDiff

torv_and_back(V::AbstractVector{<:Real}) = V, identity
torv_and_back(x::Real) = [x], V -> V[1]
torv_and_back(x::Complex) = [real(x), imag(x)], V -> Complex(V[1], V[2])

function torv_and_back(x::Tuple{Vararg{Any,N}}) where N
    xs_fs_flat = map(torv_and_back, x)
    xs = map(x -> x[1], xs_fs_flat)
    fs = map(x -> x[2], xs_fs_flat)
    offs = [0, cumsum(map(length, xs))...]
    flat_x = vcat(xs...)

    function f_back(flat_x)
        idxs = offs .+ firstindex(flat_x)
        rec_xs = getindex.(Ref(flat_x), range.(idxs[begin:end-1], idxs[begin+1:end] .- 1))
        return (map((f, x) -> f(x), fs, rec_xs)...,)
    end

    flat_x, f_back
end

function torv_and_back(x::NamedTuple{names}) where {names}
    flat_x, f_back_tpl = torv_and_back(values(x))
    f_back(flat_x) = NamedTuple{names}(f_back_tpl(flat_x))
    return flat_x, f_back
end

function torv_and_back(x::Ref)
    xval = x[]
    V, to_xval = torv_and_back(xval)
    back_to_ref(V) = Ref(to_xval(V))
    return (V, back_to_ref)
end

torv_and_back(A::AbstractArray{<:Real}) = vec(A), V -> reshape(V, size(A))

function torv_and_back(A::AbstractArray{Complex{T}, N}) where {T<:Real, N}
    RA = cat(real.(A), imag.(A), dims = N+1) 
    V, to_array = torv_and_back(RA)
    function back_to_complex(V)
        RA = to_array(V)
        Complex.(view(RA, map(_ -> :, size(A))..., 1), view(RA, map(_ -> :, size(A))..., 2))
    end
    return (V, back_to_complex)
end


function getjacobian(f, x)
    V, to_x = torv_and_back(x)
    vf(V) = torv_and_back(f(to_x(V)))[1]
    ForwardDiff.jacobian(vf, V)
end

end # !isdefined(Main, :getjacobian)

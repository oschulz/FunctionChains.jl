# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).


_zeroladj(::Type{T}) where {T<:Real} = zero(T)
_zeroladj(::Type{<:AbstractArray{T}}) where {T<:Real} = zero(T)
_zeroladj(::Type) = 0

_similar_empty(A::AbstractVector, ::Type{T}) where T = similar(A, T, 0)
_similar_empty(A::Any, ::Type{T}) where T = T[]

_push!!(A::AbstractVector{T}, x::T) where T = push!(A,x)
_push!!(A::AbstractVector{T}, x::U) where {T,U} = vcat(A, [x])

_sizehint!(A::AbstractVector, ::Union{Base.HasLength,Base.HasShape}, itr) = sizehint!(A, length(itr))
_sizehint!(A::AbstractVector, ::Base.SizeUnknown, itr) = nothing

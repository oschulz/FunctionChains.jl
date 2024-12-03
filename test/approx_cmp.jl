# This file is a part of AutoDiffOperators.jl, licensed under the MIT License (MIT).

if !isdefined(Main, :approx_cmp)

approx_cmp(a::Number, b::Number; kwargs...) = isapprox(a, b; kwargs...)
approx_cmp(a::AbstractArray{<:Number}, b::AbstractArray{<:Number}; kwargs...) = isapprox(a, b; kwargs...)
approx_cmp(a::AbstractArray, b::AbstractArray; kwargs...) = all(map((x, y) -> approx_cmp(x, y; kwargs...), a, b))
approx_cmp(a::Tuple, b::Tuple; kwargs...) = all(map((x, y) -> approx_cmp(x, y; kwargs...), a, b))
approx_cmp(a::NamedTuple{names}, b::NamedTuple{names}; kwargs...) where names = approx_cmp(values(a), values(b); kwargs...)
approx_cmp(a::NamedTuple, b::NamedTuple; kwargs...) = false

end # !isdefined(Main, :approx_cmp)

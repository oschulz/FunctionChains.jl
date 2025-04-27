# This file is a part of AutoDiffOperators.jl, licensed under the MIT License (MIT).

if !isdefined(Main, :modify_output)

modify_output(y::Number) = sign(y) * abs(y)^1.2
modify_output(y::Tuple) = map(modify_output, y)
modify_output(y::NamedTuple) = map(modify_output, y)
modify_output(y::AbstractArray) = modify_output.(y)

end # !isdefined(Main, :modify_output)

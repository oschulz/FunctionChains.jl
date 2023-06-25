# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsFlexiMapsExt

using FlexiMaps
using FunctionChains

FlexiMaps.islinear(fc::FunctionChain) = all(FlexiMaps.islinear, fc.fs)
FlexiMaps.isaffine(fc::FunctionChain) = all(FlexiMaps.isaffine, fc.fs)

end # module FunctionChainsFlexiMapsExt

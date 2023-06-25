# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsAdaptExt

using Adapt
using FunctionChains

Adapt.adapt_structure(target, f::FunctionChain) = FunctionChain(map(Base.Fix1(Adapt.adapt, target), f.fs))

end # module FunctionChainsAdaptExt

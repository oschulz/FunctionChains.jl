# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsAdaptExt

using Adapt
using FunctionChains

function Adapt.adapt_structure(target::T, fc::FunctionChain) where T
    FunctionChain(map(Base.Fix1(Adapt.adapt, target), fchainfs(fc)))
end

end # module FunctionChainsAdaptExt

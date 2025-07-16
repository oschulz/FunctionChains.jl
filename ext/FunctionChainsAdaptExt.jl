# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsAdaptExt

using Adapt

using FunctionChains
using FunctionChains: AsFunction

function Adapt.adapt_structure(target::T, fc::AsFunction) where T
    AsFunction(Adapt.adapt(target, fc.f))
end

function Adapt.adapt_structure(target::T, fc::FunctionChain) where T
    FunctionChain(map(Base.Fix1(Adapt.adapt, target), fchainfs(fc)))
end

function Adapt.adapt_structure(target::T, fp::FCartProd) where T
    FCartProd(map(Base.Fix1(Adapt.adapt, target), fcprodfs(fp)))
end

function Adapt.adapt_structure(target::T, ff::FFanout) where T
    FFanout(map(Base.Fix1(Adapt.adapt, target), ffanoutfs(ff)))
end

function Adapt.adapt_structure(target::T, fa::FAlias) where T
    falias(Adapt.adapt(target, forig(fa)))
end

end # module FunctionChainsAdaptExt

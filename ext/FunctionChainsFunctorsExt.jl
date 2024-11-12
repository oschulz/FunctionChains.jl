# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsFunctorsExt

using Functors
using FunctionChains

using FunctionChains: AsFunction

@static if !isdefined(Base, :pkgversion) || pkgversion(Functors) < v"0.5"

Functors.@functor AsFunction

Functors.@functor FunctionChain

end # Functors < v"0.5"

end # module FunctionChainsFunctorsExt

# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsFunctorsExt

using Functors
using FunctionChains

using FunctionChains: AsFunction

Functors.@functor AsFunction

Functors.@functor FunctionChain

end # module FunctionChainsFunctorsExt

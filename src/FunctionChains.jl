# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

"""
    module FunctionChains

Implements chained functions (composed functions) beyond `Base.ComposedFunction`.
"""
module FunctionChains

using ChainRulesCore
using ChangesOfVariables
using Functors
using InverseFunctions

import Adapt

include("utils.jl")
include("function_chain.jl")

end # module

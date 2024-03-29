# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

"""
    module FunctionChains

Implements chained functions (composed functions) beyond `Base.ComposedFunction`.
"""
module FunctionChains

using Tricks: static_hasmethod

include("utils.jl")
include("as_function.jl")
include("function_chain.jl")

@static if !isdefined(Base, :get_extension)
    include("../ext/FunctionChainsAdaptExt.jl")
    include("../ext/FunctionChainsChangesOfVariablesExt.jl")
    # FlexiMaps supports Julia >= v1.9 only
    include("../ext/FunctionChainsFunctorsExt.jl")
    include("../ext/FunctionChainsInverseFunctionsExt.jl")
end

end # module

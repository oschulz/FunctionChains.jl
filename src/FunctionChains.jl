# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

"""
    module FunctionChains

Implements chained functions (composed functions) beyond `Base.ComposedFunction`.
"""
module FunctionChains

using Tricks: static_hasmethod

include("utils.jl")
include("as_function.jl")
include("fchain.jl")

end # module

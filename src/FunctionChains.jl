# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

"""
    module FunctionChains

Implements chained functions (composed functions) beyond `Base.ComposedFunction`.
"""
module FunctionChains

using Base.Iterators: Take, Repeated
using Base: IteratorSize, HasShape, HasLength

using Tricks: static_hasmethod

import InverseFunctions
using InverseFunctions: inverse, NoInverse

include("utils.jl")
include("as_function.jl")
include("applyf.jl")
include("fbcast.jl")
include("fchain.jl")
include("fcomp.jl")
include("frepeat.jl")
include("fcprod.jl")

end # module

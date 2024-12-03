# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

import Test

Test.@testset "Package FunctionChains" begin
    # include("test_aqua.jl")
    include("test_applyf.jl")
    include("test_fbcast.jl")
    include("test_fchain.jl")
    include("test_fcomp.jl")
    include("test_frepeat.jl")
    include("test_fcprod.jl")
    include("test_docs.jl")
end # testset

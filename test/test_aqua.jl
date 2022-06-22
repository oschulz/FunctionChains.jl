# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import FunctionChains

Test.@testset "Aqua tests" begin
    Aqua.test_all(FunctionChains)
end # testset

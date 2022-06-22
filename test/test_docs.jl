# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

using Test
using FunctionChains
import Documenter

Documenter.DocMeta.setdocmeta!(
    FunctionChains,
    :DocTestSetup,
    :(using FunctionChains);
    recursive=true,
)
Documenter.doctest(FunctionChains)

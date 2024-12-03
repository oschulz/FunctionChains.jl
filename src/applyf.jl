# This file is a part of FunctionProducts.jl, licensed under the MIT License (MIT).


"""
    applyf(f, xs...) = f(xs...)

Function application: Apply a function `f` to one or multiple arguments.

Primarily useful to broadcast collections of function over collections of
arguments, e.g. `applyf(fs, as, bs, ...)` or to calculate derivatives in
respect to function parameters (closure contents), e.g.
`Zygote.gradient(applyf, f_withargs, xs...)`.
"""
function applyf end
export applyf

@inline applyf(f) = f()
@inline applyf(f, x) = f(x)
@inline applyf(f, a, bs...) = f(a, bs...)

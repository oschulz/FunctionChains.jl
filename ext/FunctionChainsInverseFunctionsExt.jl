# This file is a part of FunctionChains.jl, licensed under the MIT License (MIT).

module FunctionChainsInverseFunctionsExt

using InverseFunctions
using FunctionChains

_reverse(fs) = reverse(fs)
_reverse(fs::Base.Iterators.Repeated) = fs
_reverse(fs::Base.Iterators.Take{<:Base.Iterators.Repeated}) = fs

_is_noinverse(f) = f isa NoInverse
_contains_noinverse(fs) = Val(any(_is_noinverse, fs))
@generated function _contains_noinverse(fs::Tuple)
    result = any(x -> x <: NoInverse, fs.parameters)
    :(Val($result))
end
_contains_noinverse(fs::AbstractArray) = Val(false)
_contains_noinverse(fs::AbstractArray{<:NoInverse}) = Val(true)
_contains_noinverse(fs::AbstractArray{>:NoInverse}) = Val(any(_is_noinverse, fs))

function InverseFunctions.inverse(fc::FunctionChain)
    inv_fs = map(inverse, _reverse(fc.fs))
    if _contains_noinverse(inv_fs) isa Val{true}
        NoInverse(fc)
    else
        FunctionChain(inv_fs)
    end
end

InverseFunctions.inverse(f::FunctionChain{<:Base.Generator}) = NoInverse(f)


end # module FunctionChainsInverseFunctionsExt

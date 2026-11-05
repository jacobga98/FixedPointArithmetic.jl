"""
Binary fixed-point numbers for simulating fixed-point signal processing.

A fixed-point type [`Fx`](@ref) has a word length, a fraction length, a rounding mode and an overflow mode, created
with [`fixed`](@ref). Operations follow the model of hardware fixed-point types such as `ap_fixed`: products keep
full precision, sums keep the format of their operands, and rounding and overflow happen when a value is converted
to another format, which includes storing it into an array of that format.
"""
module FixedPointArithmetic

using Random: Random

export Fx,
    fixed,
    OverflowMode,
    Wrap,
    Saturate,
    SaturateSymmetric,
    rawvalue,
    fromraw,
    rawtype,
    wordlength,
    fractionlength,
    roundingmode,
    overflowmode,
    divide

include("modes.jl")
include("type.jl")
include("raw.jl")
include("conversions.jl")
include("arithmetic.jl")
include("division.jl")
include("show.jl")

end

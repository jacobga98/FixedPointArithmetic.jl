"""
    OverflowMode

How a value outside the range of a fixed-point format is brought back into it: [`Wrap`](@ref), [`Saturate`](@ref) or
[`SaturateSymmetric`](@ref). The rounding of values between two representable ones is a Base `RoundingMode`
instead, see [`fixed`](@ref).
"""
struct OverflowMode{S} end

"""
    Wrap

Overflow mode that keeps the low word-length bits of the two's complement representation, as integer hardware does
without saturation logic.
"""
const Wrap = OverflowMode{:Wrap}()

"""
    Saturate

Overflow mode that clamps to the largest or smallest representable value, `-2^(W-1)` to `2^(W-1) - 1` in units of
the least significant bit for word length `W`.
"""
const Saturate = OverflowMode{:Saturate}()

"""
    SaturateSymmetric

Overflow mode that clamps to a range symmetric around zero, `-(2^(W-1) - 1)` to `2^(W-1) - 1` in units of the least
significant bit, so that negation and `abs` never overflow. The usual choice for soft bits such as int8 LLRs.
"""
const SaturateSymmetric = OverflowMode{:SaturateSymmetric}()

# The Base rounding modes a fixed-point format can use. `RoundFromZero` is left out: it is not defined for all
# floating-point types in Base, and no fixed-point design uses it.
const ROUNDING_MODES = (
    RoundDown, RoundUp, RoundToZero, RoundNearest, RoundNearestTiesAway, RoundNearestTiesUp
)

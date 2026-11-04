module FixedPointArithmeticFixedPointNumbersExt

using FixedPointArithmetic: Fx, fromraw, rawvalue, roundingmode, Saturate
using FixedPointNumbers: Fixed

# `Fixed{T,f}` holds exactly the values of the format with the full width of T and f fraction bits.
_fxtype(::Type{Fixed{T,f}}, R) where {T,f} = Fx{8 * sizeof(T),f,R,Saturate,T}

# Into a fixed-point type, with its rounding and overflow modes.
(::Type{X})(x::Fixed{T,f}) where {X<:Fx,T,f} = X(fromraw(_fxtype(Fixed{T,f}, RoundDown), reinterpret(x)))

# Into FixedPointNumbers, with the rounding mode of `x` and saturation.
function (::Type{Q})(x::Fx) where {T,f,Q<:Fixed{T,f}}
    return reinterpret(Q, rawvalue(_fxtype(Q, roundingmode(x))(x)))
end

Base.convert(::Type{X}, x::Fixed) where {X<:Fx} = X(x)
Base.convert(::Type{Q}, x::Fx) where {Q<:Fixed} = Q(x)

end

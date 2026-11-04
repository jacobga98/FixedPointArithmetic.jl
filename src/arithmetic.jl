# Arithmetic rules:
# - Sums and differences of one format keep that format; the exact result goes through its overflow mode. Loops such
#   as `acc += x` therefore keep their type, and the accumulator width is the type the caller chose for `acc`.
# - Products are exact: the word and fraction lengths add up. Rounding happens only when the product is converted,
#   for example when it is stored into an array, so a complex product is rounded once, not per partial product.
# - Integers act as exact constants and keep the fixed-point format.
# - Floats do not mix: an operation with a float throws, so that a float constant cannot silently turn a fixed-point
#   computation into a floating-point one. Convert the constant with `X(c)` first.
# - Division and transcendental functions are not defined. Base would compute them in floating point; a fixed-point
#   algorithm has to choose its own method (a reciprocal table, CORDIC, ...), or convert explicitly with `float`.

Base.zero(::Type{X}) where {X<:Fx} = _fromraw(X, zero(rawtype(X)))
Base.zero(x::Fx) = zero(typeof(x))

# One is not representable when the format has no integer bits beside the sign, as Q15 = fixed(16, 15); the largest
# value is the closest to it.
Base.one(::Type{X}) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}} =
    _fromraw(X, _overflow(Saturate, Val(W), T, one(widen(T)) << F))
Base.one(x::Fx) = one(typeof(x))
Base.oneunit(::Type{X}) where {X<:Fx} = one(X)
Base.oneunit(x::Fx) = one(typeof(x))

Base.typemax(::Type{X}) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}} = _fromraw(X, _maxraw(Val(W), T))
Base.typemin(::Type{X}) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}} = _fromraw(X, _minraw(Val(W), T))
Base.typemin(::Type{X}) where {W,F,R,T,X<:Fx{W,F,R,SaturateSymmetric,T}} = _fromraw(X, -_maxraw(Val(W), T))
Base.typemax(x::Fx) = typemax(typeof(x))
Base.typemin(x::Fx) = typemin(typeof(x))

"""
    eps(::Type{<:Fx})

The weight `2^-F` of the least significant bit of a fixed-point format.
"""
Base.eps(::Type{X}) where {X<:Fx} = _fromraw(X, one(rawtype(X)))
Base.eps(x::Fx) = eps(typeof(x))

"""
    widen(::Type{<:Fx})

The format with twice the word length, at most 64 bits, and the same fraction length and modes.
"""
function Base.widen(::Type{Fx{W,F,R,O,T}}) where {W,F,R,O,T}
    W2 = min(2W, 64)
    return Fx{W2,F,R,O,rawtype(W2)}
end

Base.signbit(x::Fx) = signbit(x.i)
Base.iszero(x::Fx) = iszero(x.i)
Base.isfinite(::Fx) = true
Base.isnan(::Fx) = false
Base.isinf(::Fx) = false

# The value `x.i * 2^-F`, which gives Base's hashing of real numbers: `hash(fixed(16, 15)(0.5)) == hash(0.5)`.
Base.decompose(x::Fx{W,F}) where {W,F} = (x.i, -F, 1)

# Same format

Base.:+(x::X, y::X) where {X<:Fx} = _fromraw(X, _overflow(X, widen(x.i) + widen(y.i)))
Base.:-(x::X, y::X) where {X<:Fx} = _fromraw(X, _overflow(X, widen(x.i) - widen(y.i)))
Base.:-(x::X) where {X<:Fx} = _fromraw(X, _overflow(X, -widen(x.i)))
Base.abs(x::X) where {X<:Fx} = _fromraw(X, _overflow(X, abs(widen(x.i))))

Base.:(==)(x::X, y::X) where {X<:Fx} = x.i == y.i
Base.:<(x::X, y::X) where {X<:Fx} = x.i < y.i
Base.:<=(x::X, y::X) where {X<:Fx} = x.i <= y.i
Base.isless(x::X, y::X) where {X<:Fx} = x.i < y.i

# Products

function Base.:*(x::Fx{W1,F1,R,O,T1}, y::Fx{W2,F2,R,O,T2}) where {W1,F1,W2,F2,R,O,T1,T2}
    P = _producttype(Val(W1), Val(W2), Val(F1 + F2), R, O)
    T = rawtype(P)
    return _fromraw(P, (x.i % T) * (y.i % T))
end

@inline function _producttype(::Val{W1}, ::Val{W2}, ::Val{F}, R, O) where {W1,W2,F}
    W1 + W2 <= 64 || throw(
        ArgumentError(
            "the product of $W1-bit and $W2-bit fixed-point values needs $(W1 + W2) bits, more than the 64 " *
            "supported; convert an operand to a shorter format first",
        ),
    )
    return Fx{W1 + W2,F,R,O,rawtype(W1 + W2)}
end

Base.muladd(x::X, y::X, z::X) where {X<:Fx} = x * y + z
Base.muladd(x::Fx, y::Fx, z::Fx) = x * y + z

# `x^2` keeps working through `literal_pow`; other powers would change type with the exponent.
Base.literal_pow(::typeof(^), x::Fx, ::Val{2}) = x * x

# Integers as exact constants

# n·2^F in the double-width type V of the storage: exact while it is within ±2^(bits(V)-2), and clamped to that
# otherwise, which is far outside the word length. A stored value added to it cannot overflow V.
@inline function _intraw_sat(::Type{X}, n::Integer) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}}
    V = widen(T)
    m = one(V) << (_bits(V) - 2 - F)
    return (clamp(_signed(n), -m, m) % V) << F
end

# n·2^F modulo 2^bits(V), for the wrapping mode, where the arithmetic is modular anyway.
@inline _intraw_wrap(::Type{X}, n::Integer) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}} = (n % widen(T)) << F

@inline _intraw(::Type{X}, n::Integer) where {X<:Fx} = _intraw(overflowmode(X), X, n)
@inline _intraw(::OverflowMode{:Wrap}, ::Type{X}, n::Integer) where {X<:Fx} = _intraw_wrap(X, n)
@inline _intraw(::OverflowMode, ::Type{X}, n::Integer) where {X<:Fx} = _intraw_sat(X, n)

Base.:+(x::X, n::Integer) where {X<:Fx} = _fromraw(X, _overflow(X, widen(x.i) + _intraw(X, n)))
Base.:+(n::Integer, x::X) where {X<:Fx} = x + n
Base.:-(x::X, n::Integer) where {X<:Fx} = _fromraw(X, _overflow(X, widen(x.i) - _intraw(X, n)))
Base.:-(n::Integer, x::X) where {X<:Fx} = _fromraw(X, _overflow(X, _intraw(X, n) - widen(x.i)))

# x·n in the double-width type V. With saturation, an n beyond ±(2^bits(T) - 1) gives a product beyond the word
# length for any nonzero x, so clamping n there keeps the result and keeps the product inside V.
@inline function _intmul(x::X, n::Integer) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}}
    V = widen(T)
    if O === Wrap
        return _fromraw(X, _overflow(X, widen(x.i) * (n % V)))
    else
        m = (one(V) << _bits(T)) - one(V)
        return _fromraw(X, _overflow(X, widen(x.i) * (clamp(_signed(n), -m, m) % V)))
    end
end

Base.:*(x::Fx, n::Integer) = _intmul(x, n)
Base.:*(n::Integer, x::Fx) = _intmul(x, n)
Base.:*(x::Fx, b::Bool) = ifelse(b, x, zero(x))
Base.:*(b::Bool, x::Fx) = ifelse(b, x, zero(x))

for op in (:(==), :<, :<=, :isless)
    @eval begin
        Base.$op(x::X, n::Integer) where {X<:Fx} = $op(widen(x.i), _intraw_sat(X, n))
        Base.$op(n::Integer, x::X) where {X<:Fx} = $op(_intraw_sat(X, n), widen(x.i))
    end
end

# No silent floating point

@noinline function _floatmix(op, x, y)
    throw(
        ArgumentError(
            "`$op` mixes a fixed-point and a floating-point value ($(typeof(x)) and $(typeof(y))); convert the " *
            "float to the fixed-point type with `X(c)`, or the fixed-point value to float with `float(x)`",
        ),
    )
end

for op in (:+, :-, :*, :/, :(==), :<, :<=, :isless, :min, :max)
    @eval begin
        Base.$op(x::Fx, y::AbstractFloat) = _floatmix($(QuoteNode(op)), x, y)
        Base.$op(x::AbstractFloat, y::Fx) = _floatmix($(QuoteNode(op)), x, y)
    end
end

@noinline function _nomath(f, x)
    throw(
        ArgumentError(
            "`$f` is not defined for fixed-point values ($(typeof(x))): Base would compute it in floating point. " *
            "Use a fixed-point method, or convert explicitly with `float(x)`",
        ),
    )
end

for f in (
    :sqrt, :cbrt, :exp, :exp2, :exp10, :expm1, :log, :log2, :log10, :log1p, :sin, :cos, :tan, :sincos, :sinpi,
    :cospi, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh, :inv,
)
    @eval Base.$f(x::Fx) = _nomath($(QuoteNode(f)), x)
end

Base.:/(x::Fx, y::Fx) = _nomath(:/, x)
Base.:/(x::Fx, y::Integer) = _nomath(:/, x)
Base.:/(x::Integer, y::Fx) = _nomath(:/, y)
Base.atan(y::Fx, x::Fx) = _nomath(:atan, y)
Base.hypot(x::Fx, y::Fx) = _nomath(:hypot, x)
Base.:^(x::Fx, p::Integer) = _nomath(:^, x)

# Scaling by powers of two

"""
    ldexp(x::Fx, e::Integer)

`x * 2^e` in the format of `x`: a left shift goes through the overflow mode, and a right shift rounds with the
rounding mode. This is the scaling step of block floating point, for example an exponent applied to a whole block
of samples. For floating-point `x`, `ldexp` is the exact Base function, so generic code can call it for both.
"""
function Base.ldexp(x::X, e::Integer) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}}
    if e >= 0
        # Beyond W every nonzero value overflows, and the shift has to stay below the width of T.
        return _fromraw(X, _lshift(O, Val(W), x.i, Int(min(e, W))) % T)
    else
        # Beyond W + 1 the rounded result no longer depends on the shift, and the double width has room for it.
        s = Int(min(-e, W + 1))
        return _fromraw(X, _overflow(O, Val(W), T, _rshift(R, widen(x.i), s)))
    end
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{X}) where {X<:Fx}
    return _fromraw(X, rand(rng, rawvalue(typemin(X)):rawvalue(typemax(X))))
end

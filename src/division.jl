# Division and square roots, correctly rounded: the exact quotient or root is rounded once to the target format with its
# rounding mode, and brought into its range with its overflow mode. Both are integer algorithms, branch-free on the
# data and without exceptions, so they work in GPU kernels like the rest. Neither has an exact result format, unlike
# products, so `a / b` and `sqrt(a)` keep the format of `a`, and `divide` takes the format as an argument.

# The integer type with n bits for a division.
function _divtype(n::Int)
    n <= 32 && return Int32
    n <= 64 && return Int64
    n <= 128 && return Int128
    throw(
        ArgumentError(
            "the division needs $n-bit integers, more than the 128 supported; convert an operand to a shorter " *
            "format first",
        ),
    )
end

# Division by zero gives what converting ±Inf gives, and 0/0 what converting NaN gives.
@inline function _divzero(::OverflowMode{:Saturate}, ::Val{W}, ::Type{T}, n) where {W,T}
    return ifelse(n > zero(n), _maxraw(Val(W), T), ifelse(n < zero(n), _minraw(Val(W), T), zero(T)))
end
@inline function _divzero(::OverflowMode{:SaturateSymmetric}, ::Val{W}, ::Type{T}, n) where {W,T}
    return ifelse(n > zero(n), _maxraw(Val(W), T), ifelse(n < zero(n), -_maxraw(Val(W), T), zero(T)))
end
@inline _divzero(::OverflowMode{:Wrap}, ::Val{W}, ::Type{T}, n) where {W,T} = zero(T)

# The raw value of Y for the quotient (na·2^-Fa) / (nb·2^-Fb), where na fits Wa bits and nb fits Wb bits. In units of
# the result's LSB the quotient is na·2^(Fy+Fb-Fa) / nb, which one integer division rounds exactly.
@inline function _divraw(
    ::Type{Y}, na::Signed, ::Val{Wa}, ::Val{Fa}, nb::Signed, ::Val{Wb}, ::Val{Fb}
) where {Wy,Fy,Ry,Oy,Ty,Y<:Fx{Wy,Fy,Ry,Oy,Ty},Wa,Fa,Wb,Fb}
    s = Fy + Fb - Fa
    V = _divtype(max(Wa + max(s, 0), Wb + max(-s, 0)) + 1)
    n = (na % V) << max(s, 0)
    d = (nb % V) << max(-s, 0)
    dzero = iszero(d)
    q = div(n, ifelse(dzero, one(V), d), Ry)
    return ifelse(dzero, _divzero(Oy, Val(Wy), Ty, n), _overflow(Oy, Val(Wy), Ty, q))
end

# An integer operand clamped to ±2^J, in a type with room for J + 2 bits.
@inline function _clampint(n::Integer, ::Val{J}) where {J}
    I = J <= 61 ? Int64 : Int128
    m = one(I) << J
    return clamp(_signed(n), -m, m) % I
end

"""
    divide(::Type{Y}, a, b) where {Y<:Fx}

The quotient `a / b` of fixed-point values or integers in the format `Y`: the exact quotient, rounded once with the
rounding mode of `Y` and brought into its range with its overflow mode. Division by zero gives what converting ±Inf
gives (the largest or smallest value with saturation, 0 with [`Wrap`](@ref)), and 0/0 gives 0, so that `divide`
never throws.

`a / b` is `divide(typeof(a), a, b)`: the quotient keeps the format of the numerator, as sums keep the format of
their operands. With an integer numerator, `n / b` takes the format of `b`, and `inv(x)` is `1 / x`.

Integer operands far outside the range that matters for the result are clamped first, which keeps the results exact
except for quotients that overflow a [`Wrap`](@ref) format. Dividing two integers converts through `Rational` and is
meant for constants, not for kernels.

# Examples
```jldoctest
julia> Q15 = fixed(16, 15; rounding=RoundNearest);

julia> divide(Q15, Q15(0.25), Q15(0.75))
Fx{16, 15, RoundNearest, Saturate, Int16}(0.333343505859375)

julia> Q15(0.75) / Q15(0.25)  # 3 saturates in Q15
Fx{16, 15, RoundNearest, Saturate, Int16}(0.999969482421875)

julia> divide(fixed(16, 13; rounding=RoundNearest), Q15(0.75), Q15(0.25))
Fx{16, 13, RoundNearest, Saturate, Int16}(3.0)
```
"""
function divide(::Type{Y}, a::Fx{Wa,Fa}, b::Fx{Wb,Fb}) where {Y<:Fx,Wa,Fa,Wb,Fb}
    return _fromraw(Y, _divraw(Y, a.i, Val(Wa), Val(Fa), b.i, Val(Wb), Val(Fb)))
end

# An integer denominator beyond ±2^J, with J = Wa + 1 + Fy - Fa, makes the raw quotient smaller than a quarter for every
# numerator, and every rounding mode then gives what it gives at ±2^J. Clamping there keeps the result exact and the
# integers short.
function divide(::Type{Y}, a::Fx{Wa,Fa}, n::Integer) where {Wy,Fy,R,O,T,Y<:Fx{Wy,Fy,R,O,T},Wa,Fa}
    J = clamp(Wa + 1 + Fy - Fa, 1, 125)
    return _fromraw(Y, _divraw(Y, a.i, Val(Wa), Val(Fa), _clampint(n, Val(J)), Val(J + 2), Val(0)))
end

# An integer numerator beyond ±2^K, with K = (Wy - Fy) + (Wb - Fb) + 1, makes the quotient overflow the result for every
# nonzero denominator, as it does at ±2^K, so clamping keeps saturated results exact.
function divide(::Type{Y}, n::Integer, b::Fx{Wb,Fb}) where {Wy,Fy,R,O,T,Y<:Fx{Wy,Fy,R,O,T},Wb,Fb}
    K = min(Wy - Fy + Wb - Fb + 1, 125)
    return _fromraw(Y, _divraw(Y, _clampint(n, Val(K)), Val(K + 2), Val(0), b.i, Val(Wb), Val(Fb)))
end

function divide(::Type{Y}, n::Integer, m::Integer) where {W,F,R,O,T,Y<:Fx{W,F,R,O,T}}
    iszero(m) && return _fromraw(Y, _divzero(O, Val(W), T, n))
    return Y(big(n) // big(m))
end

Base.:/(a::Fx{Wa,Fa,R,O}, b::Fx{Wb,Fb,R,O}) where {Wa,Fa,Wb,Fb,R,O} = divide(typeof(a), a, b)
Base.:/(a::Fx, n::Integer) = divide(typeof(a), a, n)
Base.:/(n::Integer, b::Fx) = divide(typeof(b), n, b)
Base.inv(x::Fx) = divide(typeof(x), 1, x)

# Complex quotients as a·conj(b) / abs2(b) in the format of a: the partial products are exact and their sums keep the
# product format, so each component is rounded once. The conjugate is folded into the products, because negating the
# most negative value of b would overflow.
@inline function _cdiv(a::Complex{X}, b::Complex{<:Fx}) where {X<:Fx}
    ar, ai = reim(a)
    br, bi = reim(b)
    den = br * br + bi * bi
    return Complex(divide(X, ar * br + ai * bi, den), divide(X, ai * br - ar * bi, den))
end

function Base.:/(
    a::Complex{Fx{Wa,Fa,R,O,Ta}}, b::Complex{Fx{Wb,Fb,R,O,Tb}}
) where {Wa,Fa,Wb,Fb,R,O,Ta,Tb}
    return _cdiv(a, b)
end
# Base's method for two complex numbers of the same type would be ambiguous with the one above.
Base.:/(a::Complex{X}, b::Complex{X}) where {X<:Fx} = _cdiv(a, b)
Base.:/(a::Fx{Wa,Fa,R,O}, b::Complex{<:Fx{Wb,Fb,R,O}}) where {Wa,Fa,Wb,Fb,R,O} = _cdiv(Complex(a, zero(a)), b)

"""
    sqrt(x::Fx)

The square root of `x` in its own format, correctly rounded with its rounding mode and brought into its range with
its overflow mode (only formats without integer bits, `F == W`, can overflow). Negative values give 0, as NaN
converts to 0, so that `sqrt` never throws. The integer square root is computed bit by bit, two bits of the operand
per step, as a hardware square-root unit does.

# Examples
```jldoctest
julia> sqrt(fixed(16, 15)(0.25))
Fx{16, 15, RoundDown, Saturate, Int16}(0.5)

julia> sqrt(fixed(16, 8; rounding=RoundNearest)(2))
Fx{16, 8, RoundNearest, Saturate, Int16}(1.4140625)
```
"""
function Base.sqrt(x::X) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}}
    # In units of the result's LSB the root is sqrt(x.i·2^F), of an integer with at most W - 1 + F bits.
    B = W - 1 + F
    U = B <= 32 ? UInt32 : B <= 64 ? UInt64 : UInt128
    n = (max(x.i, zero(T)) % U) << F
    s, r = _isqrt(n, Val((B + 1) ÷ 2))
    q = s + U(_sqrtroundup(R, s, r))
    return _fromraw(X, _overflow(O, Val(W), T, q % signed(U)))
end

# floor(sqrt(n)) and the remainder n - floor(sqrt(n))^2, for n < 4^M, by the digit-by-digit method: M steps that each
# decide one bit of the root.
@inline function _isqrt(n::U, ::Val{M}) where {U<:Unsigned,M}
    root = zero(U)
    bit = one(U) << (2 * (M - 1))
    for _ in 1:M
        t = root + bit
        ge = n >= t
        n = ifelse(ge, n - t, n)
        root = ifelse(ge, (root >> 1) + bit, root >> 1)
        bit >>= 2
    end
    return root, n
end

# Whether the root rounds up from s = floor(sqrt(n)), with remainder r = n - s². The exact root is never halfway
# between two integers, so the nearest modes agree: it is at least s + 1/2 exactly when n > s² + s.
@inline _sqrtroundup(::Union{RoundingMode{:Down},RoundingMode{:ToZero}}, s, r) = false
@inline _sqrtroundup(::RoundingMode{:Up}, s, r) = !iszero(r)
@inline _sqrtroundup(
    ::Union{RoundingMode{:Nearest},RoundingMode{:NearestTiesAway},RoundingMode{:NearestTiesUp}}, s, r
) = r > s

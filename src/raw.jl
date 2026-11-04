# Integer kernels on stored values. Everything here is branch-free on the data and never throws, so that it compiles
# in GPU kernels and vectorizes on the CPU; the modes are type parameters and resolve at compile time.

# The raw range of word length W, in integer type V (which has at least W bits). For W equal to the width of V the
# shifts wrap around to typemin/typemax, which is the intended result.
@inline _maxraw(::Val{W}, ::Type{V}) where {W,V<:Signed} = (one(V) << (W - 1)) - one(V)
@inline _minraw(::Val{W}, ::Type{V}) where {W,V<:Signed} = -(one(V) << (W - 1))

_bits(::Type{V}) where {V} = 8 * sizeof(V)

# v·2^e for 0 <= e, brought into the W-bit range with overflow mode O, computed in the type V of v without widening.
# V has at least W bits.
@inline function _lshift(::OverflowMode{:Saturate}, ::Val{W}, v::V, e::Int) where {W,V<:Signed}
    hi = _maxraw(Val(W), V)
    lo = _minraw(Val(W), V)
    # Every nonzero value overflows once shifted by W or more.
    e >= W && return ifelse(v > zero(V), hi, ifelse(v < zero(V), lo, zero(V)))
    # Compare before shifting, so that the shift itself cannot overflow V.
    return ifelse(v > (hi >> e), hi, ifelse(v < (lo >> e), lo, v << e))
end

@inline function _lshift(::OverflowMode{:SaturateSymmetric}, ::Val{W}, v::V, e::Int) where {W,V<:Signed}
    hi = _maxraw(Val(W), V)
    e >= W && return ifelse(v > zero(V), hi, ifelse(v < zero(V), -hi, zero(V)))
    return ifelse(v > (hi >> e), hi, ifelse(v < -(hi >> e), -hi, v << e))
end

@inline function _lshift(::OverflowMode{:Wrap}, ::Val{W}, v::V, e::Int) where {W,V<:Signed}
    B = _bits(V)
    # Julia shifts by the width or more give 0, which is also the wrapped result of such a shift.
    w = v << e
    # Sign-extend from bit W: keep the low W bits of the two's complement representation.
    return (w << (B - W)) >> (B - W)
end

# v brought into the W-bit range with overflow mode O and narrowed to T. v can be any signed integer type, narrower or
# wider than T.
@inline function _overflow(O::OverflowMode, ::Val{W}, ::Type{T}, v::Signed) where {W,T<:Signed}
    V = promote_type(typeof(v), T)
    return _lshift(O, Val(W), v % V, 0) % T
end

@inline _overflow(::Type{<:Fx{W,F,R,O,T}}, v::Signed) where {W,F,R,O,T} = _overflow(O, Val(W), T, v)

# v / 2^s rounded with mode R, for 0 <= s < the width of V. The result is never larger in magnitude than v plus one,
# and is computed in V without overflow: floor by an arithmetic shift, then a correction from the remainder.
@inline function _rshift(R::RoundingMode, v::V, s::Int) where {V<:Signed}
    s == 0 && return v
    q = v >> s
    # The remainder of the floor division, 0 <= r < 2^s.
    r = v & ((one(V) << s) - one(V))
    half = one(V) << (s - 1)
    return q + V(_roundup(R, r, half, q, v))
end

# Whether the floor quotient q of v by 2^s, with remainder r and half = 2^(s-1), rounds up by one.
@inline _roundup(::RoundingMode{:Down}, r, half, q, v) = false
@inline _roundup(::RoundingMode{:Up}, r, half, q, v) = !iszero(r)
@inline _roundup(::RoundingMode{:ToZero}, r, half, q, v) = !iszero(r) & (v < zero(v))
@inline _roundup(::RoundingMode{:NearestTiesUp}, r, half, q, v) = r >= half
@inline _roundup(::RoundingMode{:NearestTiesAway}, r, half, q, v) = (r > half) | ((r == half) & (v >= zero(v)))
@inline _roundup(::RoundingMode{:Nearest}, r, half, q, v) = (r > half) | ((r == half) & isodd(q))

# The integer type in which a right shift by s of a value of type T is done: T itself unless the shift reaches its
# width, which `_rshift` does not handle.
_shifttype(::Type{T}, s::Int) where {T} = s < _bits(T) ? T : widen(T)

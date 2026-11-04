# Conversions into a fixed-point type apply its rounding and overflow modes and never throw, so that storing into an
# array of a fixed-point type is a quantization step, on the CPU and inside GPU kernels alike.

(::Type{X})(x::X) where {X<:Fx} = x

# Between fixed-point formats: shift by the difference in fraction bits, rounding with the target's mode when bits
# are dropped, then bring the result into the target's word length.
function (::Type{Y})(
    x::Fx{W1,F1,R1,O1,T1}
) where {W2,F2,R2,O2,T2,Y<:Fx{W2,F2,R2,O2,T2},W1,F1,R1,O1,T1}
    s = F1 - F2
    if s > 0
        V = _shifttype(T1, s)
        return _fromraw(Y, _overflow(O2, Val(W2), T2, _rshift(R2, x.i % V, s)))
    elseif s < 0
        # The shift is done in a type with room for the target word length; `_lshift` checks the range before
        # shifting, so no wider type is needed.
        V = promote_type(T1, T2)
        return _fromraw(Y, _lshift(O2, Val(W2), x.i % V, -s) % T2)
    else
        return _fromraw(Y, _overflow(O2, Val(W2), T2, x.i))
    end
end

"""
    (::Type{X})(x::AbstractFloat) where {X<:Fx}

Quantize `x` to the fixed-point type `X` with its rounding and overflow modes. The rounding is exact with respect to
the value of `x`. NaN becomes 0, and so do ±Inf with [`Wrap`](@ref).
"""
function (::Type{X})(x::AbstractFloat) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}}
    v = round(ldexp(_workfloat(x), F), R)
    return _fromraw(X, _floatraw(O, Val(W), T, v))
end

_workfloat(x::Union{Float32,Float64,BigFloat}) = x
# Exact for Float16 and BFloat16, whose range is too small for the scaling by 2^F.
_workfloat(x::AbstractFloat) = Float32(x)

_ftrunc(::Type{I}, c::Base.IEEEFloat) where {I} = unsafe_trunc(I, c)
_ftrunc(::Type{I}, c::AbstractFloat) where {I} = trunc(I, c)

# The raw value of word length W for the integer-valued float v.
@inline function _floatraw(O::OverflowMode, ::Val{W}, ::Type{T}, v::AbstractFloat) where {W,T}
    # Clamping to ±2^(W-1), which is exact in floating point and fits a W+1 bit integer, leaves the exact saturation
    # to the integer step.
    lim = ldexp(one(v), W - 1)
    c = ifelse(isnan(v), zero(v), clamp(v, -lim, lim))
    return _overflow(O, Val(W), T, _ftrunc(W < 64 ? rawtype(W + 1) : Int128, c))
end

@inline function _floatraw(O::OverflowMode{:Wrap}, ::Val{W}, ::Type{T}, v::AbstractFloat) where {W,T}
    # Wrapping needs the low bits, so the value goes through a wide integer. Values beyond its range, which only
    # occur for absurd inputs, saturate before wrapping.
    I = W <= 32 ? Int64 : Int128
    lim = ldexp(one(v), _bits(I) - 2)
    c = ifelse(isfinite(v), clamp(v, -lim, lim), zero(v))
    return _overflow(O, Val(W), T, _ftrunc(I, c))
end

# Integers: the value is first brought into the storage type, saturating or wrapping, which leaves the result of the
# shift by F unchanged: a value that does not fit the storage type overflows the word length anyway.
function (::Type{X})(n::Integer) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}}
    return _fromraw(X, _lshift(O, Val(W), _tostorage(O, T, n), F))
end

_tostorage(::OverflowMode{:Wrap}, ::Type{T}, n::Integer) where {T} = n % T
_tostorage(::OverflowMode, ::Type{T}, n::Integer) where {T} = clamp(_signed(n), typemin(T), typemax(T)) % T

# Exact constants, mostly for coefficients computed at construction time; not meant for kernels.
function (::Type{X})(x::Rational) where {W,F,R,O,T,X<:Fx{W,F,R,O,T}}
    v = round(BigInt, x * big(2)^F, R)
    O === Wrap && return _fromraw(X, _overflow(O, Val(W), T, v % T))
    return _fromraw(X, _overflow(O, Val(W), T, Int128(clamp(v, typemin(Int128), typemax(Int128)))))
end

(::Type{X})(x::AbstractIrrational) where {X<:Fx} = X(BigFloat(x; precision=256))

# To floating point. The conversion is exact when the float type has at least W significand bits.
(::Type{T})(x::Fx{W,F}) where {T<:AbstractFloat,W,F} = ldexp(T(x.i), -F)

Base.float(::Type{<:Fx{W}}) where {W} = W <= 24 ? Float32 : Float64
Base.float(x::Fx) = float(typeof(x))(x)
Base.AbstractFloat(x::Fx) = float(x)

# Mixed formats with the same modes promote to the format that holds both exactly: the most fraction bits and the
# most integer bits. Formats with different modes do not promote, so mixing them needs an explicit conversion.
function Base.promote_rule(
    ::Type{Fx{W1,F1,R,O,T1}}, ::Type{Fx{W2,F2,R,O,T2}}
) where {W1,F1,W2,F2,R,O,T1,T2}
    F = max(F1, F2)
    W = max(W1 - F1, W2 - F2) + F
    return W <= 64 ? Fx{W,F,R,O,rawtype(W)} : Union{}
end

# Integers convert into the fixed-point format, for generic code such as `clamp(x, -1, 1)`. Arithmetic and
# comparisons with integers have their own exact methods (see arithmetic.jl) and do not go through this.
Base.promote_rule(::Type{X}, ::Type{<:Integer}) where {X<:Fx} = X

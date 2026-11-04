"""
    Fx{W,F,R,O,T} <: Real

A signed binary fixed-point number with word length `W` bits and `F` fraction bits: the value `i / 2^F` of an integer
`i` in the `W`-bit two's complement range. `R` is the rounding mode (a Base `RoundingMode`) and `O` the
[`OverflowMode`](@ref) applied when a value is converted to this format. `T` is the smallest signed integer type
with at least `W` bits, which stores `i`. Create the type with [`fixed`](@ref) rather than spelling out the
parameters.

`W` is at most 64, and `F` is between 0 and `W`.

# Examples
```jldoctest
julia> Q15 = fixed(16, 15)
Fx{16, 15, RoundDown, Saturate, Int16}

julia> x = Q15(0.3)
Fx{16, 15, RoundDown, Saturate, Int16}(0.29998779296875)

julia> rawvalue(x)
9830
```
"""
struct Fx{W,F,R,O,T<:Signed} <: Real
    i::T

    # The only constructor that takes a raw integer; `fromraw` is the public way in. `X(n::Integer)` converts the
    # value `n` instead.
    Fx{W,F,R,O,T}(i::T, ::Val{:raw}) where {W,F,R,O,T<:Signed} = new{W,F,R,O,T}(i)
end

# `i` must already be of the storage type and within the word length.
@inline _fromraw(::Type{X}, i::Signed) where {X<:Fx} = X(i, Val(:raw))

"""
    rawtype(W::Integer)
    rawtype(::Type{<:Fx})

The signed integer type that stores a word length of `W` bits: the smallest of `Int8`, `Int16`, `Int32` and `Int64`
that holds `W` bits.
"""
function rawtype(W::Integer)
    W <= 8 && return Int8
    W <= 16 && return Int16
    W <= 32 && return Int32
    W <= 64 && return Int64
    throw(ArgumentError("word lengths above 64 bits are not supported, got $W"))
end

"""
    fixed(W::Integer, F::Integer; rounding=RoundDown, overflow=Saturate)

The fixed-point type [`Fx`](@ref) with word length `W` bits, of which `F` are fraction bits, so that its least
significant bit has weight `2^-F`. `rounding` is one of `RoundDown` (floor, what an arithmetic right shift does),
`RoundUp`, `RoundToZero`, `RoundNearest` (ties to even), `RoundNearestTiesAway` and `RoundNearestTiesUp` (round
half up, what adding half an LSB before the shift does). `overflow` is [`Wrap`](@ref), [`Saturate`](@ref) or
[`SaturateSymmetric`](@ref).

Both modes apply when a value is converted into this type, including when it is stored into an array of it.

# Examples
```jldoctest
julia> fixed(8, 0; overflow=SaturateSymmetric)  # int8 bits in -127:127
Fx{8, 0, RoundDown, SaturateSymmetric, Int8}

julia> fixed(6, 2; rounding=RoundNearest)        # 6-bit values with step 0.25
Fx{6, 2, RoundNearest, Saturate, Int8}
```
"""
function fixed(W::Integer, F::Integer; rounding::RoundingMode=RoundDown, overflow::OverflowMode=Saturate)
    1 <= W <= 64 || throw(ArgumentError("the word length must be between 1 and 64 bits, got $W"))
    0 <= F <= W || throw(ArgumentError("the fraction length must be between 0 and the word length $W, got $F"))
    rounding in ROUNDING_MODES ||
        throw(ArgumentError("rounding mode $rounding is not supported, use one of $ROUNDING_MODES"))
    return Fx{Int(W),Int(F),rounding,overflow,rawtype(W)}
end

"""
    wordlength(::Type{<:Fx}), wordlength(x::Fx)

The word length in bits of a fixed-point type or value.
"""
wordlength(::Type{<:Fx{W}}) where {W} = W

"""
    fractionlength(::Type{<:Fx}), fractionlength(x::Fx)

The number of fraction bits of a fixed-point type or value: its least significant bit has weight `2^-F`.
"""
fractionlength(::Type{<:Fx{W,F}}) where {W,F} = F

"""
    roundingmode(::Type{<:Fx}), roundingmode(x::Fx)

The rounding mode applied when a value is converted to this fixed-point type.
"""
roundingmode(::Type{<:Fx{W,F,R}}) where {W,F,R} = R

"""
    overflowmode(::Type{<:Fx}), overflowmode(x::Fx)

The [`OverflowMode`](@ref) applied when a value is converted to this fixed-point type.
"""
overflowmode(::Type{<:Fx{W,F,R,O}}) where {W,F,R,O} = O

rawtype(::Type{<:Fx{W,F,R,O,T}}) where {W,F,R,O,T} = T

for f in (:wordlength, :fractionlength, :roundingmode, :overflowmode, :rawtype)
    @eval $f(x::Fx) = $f(typeof(x))
end

"""
    rawvalue(x::Fx)

The stored integer `i` of `x`, whose value is `i / 2^F`.
"""
rawvalue(x::Fx) = x.i

"""
    fromraw(::Type{X}, i::Integer) where {X<:Fx}

The value of fixed-point type `X` with stored integer `i`, for example to read integer samples from a file or from
another implementation. An `i` outside the word length is brought into it with the overflow mode of `X`.

# Examples
```jldoctest
julia> fromraw(fixed(16, 15), 16384)
Fx{16, 15, RoundDown, Saturate, Int16}(0.5)
```
"""
fromraw(::Type{X}, i::Integer) where {X<:Fx} = _fromraw(X, _overflow(X, _signed(i)))

_signed(i::Signed) = i
_signed(i::Unsigned) = _signed_unsigned(i)
_signed(i::Bool) = Int8(i)
# An unsigned integer whose top bit is set does not fit the signed type of the same width.
_signed_unsigned(i::UInt8) = Int16(i)
_signed_unsigned(i::UInt16) = Int32(i)
_signed_unsigned(i::UInt32) = Int64(i)
_signed_unsigned(i::Union{UInt64,UInt128}) = Int128(min(i, typemax(Int128)))

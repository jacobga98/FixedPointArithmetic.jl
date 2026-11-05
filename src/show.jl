Base.show(io::IO, ::OverflowMode{S}) where {S} = print(io, S)

# The Base rounding modes print as `RoundingMode{:Down}()`; inside the type they print as the constant names, so that
# the printed type is valid Julia that evaluates to it.
_roundname(::RoundingMode{S}) where {S} = "Round" * String(S)

function Base.show(io::IO, X::Type{<:Fx})
    # Types with free parameters, such as those in method signatures, print the Base way.
    if X isa DataType && roundingmode(X) isa RoundingMode
        print(io, "Fx{", wordlength(X), ", ", fractionlength(X), ", ", _roundname(roundingmode(X)), ", ")
        print(io, overflowmode(X), ", ", rawtype(X), "}")
    else
        invoke(show, Tuple{IO,Type}, io, X)
    end
    return nothing
end

function Base.show(io::IO, x::X) where {X<:Fx}
    # Float64 holds every value exactly up to 53-bit words.
    if get(io, :compact, false) || get(io, :typeinfo, Any) === X
        print(io, Float64(x))
    else
        print(io, X, "(", Float64(x), ")")
    end
    return nothing
end

function Base.show(io::IO, z::Complex{X}) where {X<:Fx}
    w = Complex{Float64}(z)
    if get(io, :compact, false) || get(io, :typeinfo, Any) === Complex{X}
        show(io, w)
    else
        print(io, Complex{X}, "(", w, ")")
    end
    return nothing
end

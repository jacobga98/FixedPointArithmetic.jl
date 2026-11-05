using FixedPointArithmetic
using Random
using Test
import FixedPointNumbers as FPN

const ROUNDINGS = (RoundDown, RoundUp, RoundToZero, RoundNearest, RoundNearestTiesAway, RoundNearestTiesUp)
const OVERFLOWS = (Wrap, Saturate, SaturateSymmetric)

# Reference model in exact arithmetic

exact(x::Fx) = rawvalue(x) // big(2)^fractionlength(x)

function refoverflow(X, v::Integer)
    W = wordlength(X)
    hi = big(2)^(W - 1) - 1
    O = overflowmode(X)
    O === Saturate && return clamp(big(v), -hi - 1, hi)
    O === SaturateSymmetric && return clamp(big(v), -hi, hi)
    return mod(big(v) + hi + 1, big(2)^W) - hi - 1
end

refraw(X, r::Real) = refoverflow(X, round(BigInt, r * big(2)^fractionlength(X), roundingmode(X)))

scale2(r, e) = e >= 0 ? r * big(2)^e : r / big(2)^(-e)

allvalues(X) = [fromraw(X, i) for i in rawvalue(typemin(X)):rawvalue(typemax(X))]
samples(rng, X, n) = wordlength(X) <= 8 ? allvalues(X) : vcat(typemin(X), typemax(X), zero(X), rand(rng, X, n))

# Up to five items for which `pred` fails, so that a failing test shows examples.
failures(pred, items) = collect(Iterators.take(Iterators.filter(!pred, items), 5))

const RNG = Xoshiro(20261101)

@testset "FixedPointArithmetic" begin
    @testset "formats" begin
        @test fixed(16, 15) === Fx{16,15,RoundDown,Saturate,Int16}
        @test fixed(6, 2; rounding=RoundNearest, overflow=Wrap) === Fx{6,2,RoundNearest,Wrap,Int8}
        @test rawtype(fixed(6, 2)) === Int8
        @test rawtype(fixed(17, 0)) === Int32
        @test rawtype(fixed(64, 0)) === Int64
        @test_throws ArgumentError fixed(0, 0)
        @test_throws ArgumentError fixed(65, 0)
        @test_throws ArgumentError fixed(8, 9)
        @test_throws ArgumentError fixed(8, -1)
        @test_throws ArgumentError fixed(8, 4; rounding=RoundFromZero)
        @test isbitstype(fixed(16, 15))
        @test sizeof(fixed(6, 2)) == 1
        X = fixed(12, 3; rounding=RoundNearest, overflow=Wrap)
        @test (wordlength(X), fractionlength(X), roundingmode(X), overflowmode(X)) == (12, 3, RoundNearest, Wrap)
        @test wordlength(X(1)) == 12
    end

    @testset "printing" begin
        for X in (fixed(16, 15), fixed(6, 2; rounding=RoundNearestTiesUp, overflow=SaturateSymmetric))
            @test eval(Meta.parse(string(X))) === X
        end
        @test repr(Fx) == "Fx"
        @test repr(Fx{16,15}) == "Fx{16, 15}"
        @test repr(Vector{<:Fx}) == "Vector{<:Fx}"
        @test sprint(show, Tuple{Fx{W,F,R,O,T}} where {W,F,R,O,T}) isa String
        Q15 = fixed(16, 15)
        @test repr(Q15(0.5)) == "Fx{16, 15, RoundDown, Saturate, Int16}(0.5)"
        @test repr([Q15(0.5), Q15(-0.25)]) == "Fx{16, 15, RoundDown, Saturate, Int16}[0.5, -0.25]"
        @test repr(Complex(Q15(0.5), Q15(-0.25))) == "Complex{Fx{16, 15, RoundDown, Saturate, Int16}}(0.5 - 0.25im)"
    end

    @testset "from floats" begin
        for (W, F) in ((6, 2), (8, 0), (16, 15), (32, 30)), R in ROUNDINGS, O in OVERFLOWS
            X = fixed(W, F; rounding=R, overflow=O)
            lsb = 2.0^-F
            span = 2.0^(W - 1 - F)
            # Random values over twice the range, exact ties, and values at the range edges
            xs = vcat(
                (rand(RNG, 300) .- 0.5) .* 4span,
                (rand(RNG, -2^(W - 1):(2^(W - 1) - 1), 50) .+ 0.5) .* lsb,
                [0.0, -0.0, lsb / 2, -lsb / 2, span, -span, span - lsb, -span - lsb, 3span, -3span],
            )
            @test isempty(failures(x -> rawvalue(X(x)) == refraw(X, Rational{BigInt}(x)), xs))
            xs32 = Float32.(xs)
            @test isempty(failures(x -> rawvalue(X(x)) == refraw(X, Rational{BigInt}(x)), xs32))
            @test rawvalue(X(Float16(0.3))) == refraw(X, Rational{BigInt}(Float16(0.3)))
            @test rawvalue(X(big"0.3")) == refraw(X, Rational{BigInt}(big"0.3"))

            @test iszero(X(NaN))
            @test X(Inf) == (O === Wrap ? zero(X) : typemax(X))
            @test X(-Inf) == (O === Wrap ? zero(X) : typemin(X))
        end
    end

    @testset "from integers, rationals and irrationals" begin
        ns = vcat(-300:300, [typemin(Int64), typemax(Int64), big(10)^30, -big(10)^30])
        others = (true, false, UInt8(200), typemax(UInt64), Int8(-128))
        for (W, F) in ((6, 2), (8, 0), (16, 15), (16, 4), (32, 8)), O in OVERFLOWS
            X = fixed(W, F; overflow=O)
            @test isempty(failures(n -> rawvalue(X(n)) == refoverflow(X, big(n) * big(2)^F), ns))
            @test isempty(failures(n -> rawvalue(X(n)) == refoverflow(X, big(n) * big(2)^F), others))
            for R in ROUNDINGS
                Y = fixed(W, F; rounding=R, overflow=O)
                for r in (1 // 3, -7 // 5, 22 // 7, -1 // 2^(F + 1))
                    @test rawvalue(Y(r)) == refraw(Y, r)
                end
                @test rawvalue(Y(π)) == refraw(Y, Rational{BigInt}(BigFloat(π; precision=256)))
            end
        end
    end

    @testset "between formats" begin
        sources = (fixed(8, 7; overflow=Wrap), fixed(8, 0), fixed(16, 3), fixed(32, 24))
        targets = ((6, 2), (8, 4), (4, 0), (8, 8), (16, 15), (32, 20))
        for S in sources, (W, F) in targets, R in ROUNDINGS, O in OVERFLOWS
            Y = fixed(W, F; rounding=R, overflow=O)
            xs = samples(RNG, S, 300)
            @test isempty(failures(x -> rawvalue(Y(x)) == refraw(Y, exact(x)), xs))
        end
        # Storing into an array converts
        Y = fixed(8, 4; rounding=RoundNearest)
        xs = samples(RNG, fixed(16, 10), 100)
        ys = Vector{Y}(undef, length(xs))
        ys .= xs
        @test ys == Y.(xs)
        @test convert(Y, xs[5]) === Y(xs[5])
    end

    @testset "sums" begin
        for (W, F) in ((6, 3), (8, 0), (16, 15), (32, 31)), O in OVERFLOWS
            X = fixed(W, F; overflow=O)
            xs = samples(RNG, X, 200)
            pairs = W <= 8 ? vec(collect(Iterators.product(xs, xs))) : collect(zip(xs, reverse(xs)))
            @test isempty(failures(((a, b),) -> rawvalue(a + b) == refoverflow(X, big(a.i) + b.i), pairs))
            @test isempty(failures(((a, b),) -> rawvalue(a - b) == refoverflow(X, big(a.i) - b.i), pairs))
            @test isempty(failures(a -> rawvalue(-a) == refoverflow(X, -big(a.i)), xs))
            @test isempty(failures(a -> rawvalue(abs(a)) == refoverflow(X, abs(big(a.i))), xs))
            @test typeof(xs[1] + xs[2]) === X
        end
        X = fixed(8, 0; overflow=SaturateSymmetric)
        @test rawvalue(fromraw(X, -128)) == -127
        @test rawvalue(abs(typemin(X))) == 127
    end

    @testset "selection" begin
        X = fixed(6, 3)
        pairs = vec(collect(Iterators.product(allvalues(X), allvalues(X))))
        @test isempty(failures(((a, b),) -> min(a, b) === (exact(a) <= exact(b) ? a : b), pairs))
        @test isempty(failures(((a, b),) -> max(a, b) === (exact(a) >= exact(b) ? a : b), pairs))
        @test ifelse(true, X(1), X(2)) === X(1)
        @test ifelse(false, X(1), X(2)) === X(2)
        @test clamp(X(3), X(-1), X(1)) === X(1)
    end

    @testset "products" begin
        for R in (RoundDown, RoundNearest), O in OVERFLOWS
            X1 = fixed(6, 3; rounding=R, overflow=O)
            X2 = fixed(5, 1; rounding=R, overflow=O)
            P = fixed(11, 4; rounding=R, overflow=O)
            pairs = vec(collect(Iterators.product(allvalues(X1), allvalues(X2))))
            @test all(((a, b),) -> typeof(a * b) === P, pairs)
            @test isempty(failures(((a, b),) -> exact(a * b) == exact(a) * exact(b), pairs))
            for W in (16, 32)
                X = fixed(W, W - 1; rounding=R, overflow=O)
                xs = samples(RNG, X, 500)
                @test isempty(failures(((a, b),) -> exact(a * b) == exact(a) * exact(b), zip(xs, reverse(xs))))
            end
        end
        @test_throws ArgumentError fixed(40, 0)(1) * fixed(40, 0)(1)
        @test fixed(16, 15)(0.5)^2 === fixed(16, 15)(0.5) * fixed(16, 15)(0.5)
        @test muladd(fixed(8, 4)(1), fixed(8, 4)(2), fixed(8, 4)(3)) === fixed(16, 8)(5)
    end

    @testset "complex products round once" begin
        for R in ROUNDINGS, O in OVERFLOWS
            X = fixed(16, 15; rounding=R, overflow=O)
            P = fixed(32, 30; rounding=R, overflow=O)
            m1 = fromraw(X, -32768)
            cases = [(Complex(m1, m1), Complex(m1, -m1)); collect(zip(rand(RNG, Complex{X}, 300), rand(RNG, Complex{X}, 300)))]
            for (a, b) in cases
                z = a * b
                @test typeof(z) === Complex{P}
                # The partial products are exact; their sum keeps the 32-bit format
                re = refoverflow(P, big(real(a).i) * real(b).i - big(imag(a).i) * imag(b).i)
                im_ = refoverflow(P, big(real(a).i) * imag(b).i + big(imag(a).i) * real(b).i)
                @test rawvalue(real(z)) == re && rawvalue(imag(z)) == im_
                # and rounding happens once, when the product is stored in the 16-bit format
                c = Complex{X}(z)
                @test rawvalue(real(c)) == refraw(X, re // big(2)^30)
                @test rawvalue(imag(c)) == refraw(X, im_ // big(2)^30)
            end
            a = rand(RNG, Complex{X})
            @test abs2(a) === real(a) * real(a) + imag(a) * imag(a)
        end
    end

    @testset "mixed formats" begin
        P = fixed(19, 15)
        @test promote_type(fixed(8, 4), fixed(16, 15)) === P
        for (a, b) in zip(samples(RNG, fixed(8, 4), 0), rand(RNG, fixed(16, 15), 256))
            # The promoted format holds both operands exactly; the sum keeps it, so it can saturate
            @test typeof(a + b) === P
            @test rawvalue(a + b) == refoverflow(P, numerator((exact(a) + exact(b)) * big(2)^15))
            @test (a < b) == (exact(a) < exact(b))
        end
        @test_throws ErrorException fixed(8, 4)(1) + fixed(8, 4; overflow=Wrap)(1)
    end

    @testset "integer constants" begin
        ns = (-1000, -3, -1, 0, 1, 2, 3, 200, typemax(Int64), typemin(Int64), big(10)^25, true, false, UInt8(255))
        for (W, F) in ((8, 3), (16, 15), (16, 0), (32, 8)), O in OVERFLOWS
            X = fixed(W, F; overflow=O)
            xs = samples(RNG, X, 30)
            k = big(2)^F
            for n in ns
                @test isempty(failures(x -> rawvalue(x + n) == refoverflow(X, x.i + n * k), xs))
                @test isempty(failures(x -> rawvalue(n + x) == refoverflow(X, x.i + n * k), xs))
                @test isempty(failures(x -> rawvalue(x - n) == refoverflow(X, x.i - n * k), xs))
                @test isempty(failures(x -> rawvalue(n - x) == refoverflow(X, n * k - x.i), xs))
                @test isempty(failures(x -> rawvalue(x * n) == refoverflow(X, big(x.i) * n), xs))
                @test isempty(failures(x -> rawvalue(n * x) == refoverflow(X, big(x.i) * n), xs))
                for op in (<, <=, ==, isless)
                    @test isempty(failures(x -> op(x, n) == op(exact(x), n), xs))
                    @test isempty(failures(x -> op(n, x) == op(n, exact(x)), xs))
                end
            end
        end
        Q15 = fixed(16, 15)
        @test clamp(Q15(0.7), -1, 0) === zero(Q15)
    end

    @testset "ldexp" begin
        for R in ROUNDINGS, O in OVERFLOWS
            X = fixed(16, 8; rounding=R, overflow=O)
            xs = samples(RNG, X, 200)
            for e in -20:20
                @test isempty(failures(x -> rawvalue(ldexp(x, e)) == refraw(X, scale2(exact(x), e)), xs))
            end
            @test ldexp(X(3), typemax(Int)) == (O === Wrap ? zero(X) : typemax(X))
        end
    end

    @testset "constants" begin
        Q15 = fixed(16, 15)
        @test iszero(zero(Q15))
        @test one(Q15) === typemax(Q15)
        @test exact(one(fixed(8, 4))) == 1
        @test exact(eps(Q15)) == 1 // 2^15
        @test rawvalue(typemin(fixed(8, 0))) == -128
        @test rawvalue(typemin(fixed(8, 0; overflow=SaturateSymmetric))) == -127
        @test rawvalue(typemax(fixed(6, 0))) == 31
        @test widen(Q15) === fixed(32, 15)
        @test float(Q15) === Float32
        @test float(fixed(32, 31)) === Float64
        @test Float64(Q15(0.25)) === 0.25
        @test sign(Q15(-0.5)) === Q15(-1)
        @test signbit(Q15(-0.5)) && !signbit(Q15(0.5))
    end

    @testset "division" begin
        # The exact quotient rounded once; division by zero as ±Inf, and 0/0 as NaN, convert
        refdiv(Y, p, q) = iszero(q) ? (iszero(p) ? big(0) : big(rawvalue(Y(sign(p) * Inf)))) : refraw(Y, p // q)

        for R in ROUNDINGS, O in OVERFLOWS
            A = fixed(6, 3; rounding=R, overflow=O)
            B = fixed(5, 2; rounding=R, overflow=O)
            pairs = vec(collect(Iterators.product(allvalues(A), allvalues(B))))
            for Y in (A, fixed(8, 6; rounding=R, overflow=O), fixed(4, 0; rounding=R, overflow=O))
                @test isempty(
                    failures(((a, b),) -> rawvalue(divide(Y, a, b)) == refdiv(Y, exact(a), exact(b)), pairs)
                )
            end
            @test all(((a, b),) -> a / b === divide(A, a, b), pairs)
        end

        for R in (RoundDown, RoundNearest, RoundNearestTiesUp), O in OVERFLOWS
            f(W, F) = fixed(W, F; rounding=R, overflow=O)
            formats = (
                (f(16, 15), f(16, 15), f(16, 15)),
                (f(16, 15), f(16, 4), f(16, 10)),
                (f(32, 30), f(32, 30), f(16, 15)),
                (f(32, 8), f(16, 15), f(32, 16)),
                (f(64, 32), f(32, 16), f(64, 40)),  # needs 128-bit integers
            )
            for (A, B, Y) in formats
                pairs = zip(samples(RNG, A, 300), samples(RNG, B, 300))
                @test isempty(
                    failures(((a, b),) -> rawvalue(divide(Y, a, b)) == refdiv(Y, exact(a), exact(b)), pairs)
                )
            end
        end

        ns = (-1000, -3, -1, 0, 1, 2, 3, 7, typemax(Int64), typemin(Int64), big(10)^25, true, false, UInt8(255))
        for (W, F) in ((8, 3), (16, 15), (16, 0), (32, 8)), R in (RoundDown, RoundNearest), O in OVERFLOWS
            X = fixed(W, F; rounding=R, overflow=O)
            xs = samples(RNG, X, 30)
            for n in ns
                @test isempty(failures(x -> rawvalue(x / n) == refdiv(X, exact(x), big(n) // 1), xs))
                # Large integer numerators are clamped, which changes only quotients that overflow with Wrap
                if O !== Wrap || -3 <= n <= 3
                    @test isempty(failures(x -> rawvalue(n / x) == refdiv(X, big(n) // 1, exact(x)), xs))
                end
            end
            @test inv(xs[5]) === 1 / xs[5]
        end

        Y = fixed(16, 13; rounding=RoundNearest)
        @test rawvalue(divide(Y, 1, 3)) == refraw(Y, 1 // 3)
        @test divide(Y, 5, 0) === typemax(Y)
        @test divide(Y, -5, 0) === typemin(Y)
        @test iszero(divide(Y, 0, 0))
        for O in OVERFLOWS
            X = fixed(16, 15; overflow=O)
            @test X(0.5) / zero(X) === (O === Wrap ? zero(X) : typemax(X))
            @test X(-0.5) / 0 === (O === Wrap ? zero(X) : typemin(X))
            @test iszero(zero(X) / zero(X))
        end
        @test_throws ErrorException fixed(8, 4)(1) / fixed(8, 4; overflow=Wrap)(1)
    end

    @testset "complex division rounds once" begin
        refdiv(Y, p, q) = iszero(q) ? (iszero(p) ? big(0) : big(rawvalue(Y(sign(p) * Inf)))) : refraw(Y, p // q)
        for R in (RoundDown, RoundNearest), O in OVERFLOWS
            X = fixed(16, 15; rounding=R, overflow=O)
            P = fixed(32, 30; rounding=R, overflow=O)
            m1 = fromraw(X, -32768)
            cases = [
                (Complex(m1, m1), Complex(m1, m1))
                (Complex(X(0.5), X(0.25)), zero(Complex{X}))
                collect(zip(rand(RNG, Complex{X}, 300), rand(RNG, Complex{X}, 300)))
            ]
            for (a, b) in cases
                c = a / b
                @test typeof(c) === Complex{X}
                ar, ai, br, bi = big.(rawvalue.((real(a), imag(a), real(b), imag(b))))
                # a·conj(b) and abs2(b) are exact sums of products in the 32-bit product format
                re = refoverflow(P, ar * br + ai * bi)
                im_ = refoverflow(P, ai * br - ar * bi)
                den = refoverflow(P, br^2 + bi^2)
                @test rawvalue(real(c)) == refdiv(X, re // 1, den // 1)
                @test rawvalue(imag(c)) == refdiv(X, im_ // 1, den // 1)
            end
            a = rand(RNG, Complex{X})
            x = rand(RNG, X)
            @test x / a === Complex(x, zero(x)) / a
            @test a / x === Complex(real(a) / x, imag(a) / x)
            @test a / 2 === Complex(real(a) / 2, imag(a) / 2)
        end
    end

    @testset "sqrt" begin
        # floor(sqrt(N)) for N = x.i·2^F, rounded up when the root is at least s + 1/2, that is when 4N > (2s + 1)²
        function refsqrt(X, x)
            rawvalue(x) < 0 && return big(0)
            N = big(rawvalue(x)) << fractionlength(X)
            s = isqrt(N)
            R = roundingmode(X)
            up = R === RoundDown || R === RoundToZero ? false : R === RoundUp ? s^2 != N : 4N > (2s + 1)^2
            return refoverflow(X, s + up)
        end
        for (W, F) in ((6, 2), (8, 0), (8, 4), (8, 7), (8, 8)), R in ROUNDINGS, O in OVERFLOWS
            X = fixed(W, F; rounding=R, overflow=O)
            @test isempty(failures(x -> rawvalue(sqrt(x)) == refsqrt(X, x), allvalues(X)))
        end
        for (W, F) in ((16, 15), (16, 4), (32, 30), (32, 0), (64, 32), (64, 63)), R in ROUNDINGS
            X = fixed(W, F; rounding=R)
            @test isempty(failures(x -> rawvalue(sqrt(x)) == refsqrt(X, x), samples(RNG, X, 300)))
        end
        @test sqrt(fixed(16, 15)(0.25)) === fixed(16, 15)(0.5)
        @test iszero(sqrt(fixed(16, 15)(-0.25)))
        @test abs(Float64(sqrt(fixed(32, 16; rounding=RoundNearest)(2))) - sqrt(2)) <= 2.0^-17
    end

    @testset "no silent floating point" begin
        x = fixed(16, 15)(0.5)
        for f in (y -> y + 0.5, y -> 0.5 + y, y -> y * 0.5f0, y -> y / 0.5, y -> y < 0.1, y -> y == 0.5, y -> max(y, 0.1))
            @test_throws ArgumentError f(x)
        end
        for f in (exp, log, sin, cos, atan, y -> y^3)
            @test_throws ArgumentError f(x)
        end
        @test_throws ArgumentError cis(x)
        # Base's complex functions would convert to float first
        z = Complex(x, x)
        for f in (sqrt, exp, log, tan, atan, asin, exp2, log10, abs, angle, w -> w^3)
            @test_throws ArgumentError f(z)
        end
        @test z^2 === z * z
    end

    @testset "arrays" begin
        Q15 = fixed(16, 15)
        a = rand(RNG, Q15, 1000)
        b = rand(RNG, Q15, 1000)
        c = similar(a)
        c .= a .* b
        @test c == [Q15(x * y) for (x, y) in zip(a, b)]
        # Saturating sums depend on the order; below Base's pairwise block size `sum` adds left to right
        @test sum(a[1:100]) === foldl(+, a[1:100])
        @test typeof(sum(a)) === Q15
        r = rand(RNG, fixed(6, 2; overflow=SaturateSymmetric), 1000)
        @test all(x -> -31 <= rawvalue(x) <= 31, r)
    end

    @testset "hashing" begin
        @test hash(fixed(16, 15)(0.5)) == hash(0.5) == hash(fixed(8, 4)(0.5))
        @test isequal(fixed(16, 15)(0.5), fixed(8, 4)(0.5))
        @test length(Set([fixed(16, 15)(0.5), fixed(8, 4)(0.5)])) == 1
    end

    @testset "FixedPointNumbers" begin
        q = FPN.Q0f15(0.3)
        @test rawvalue(fixed(16, 15)(q)) == reinterpret(q)
        @test rawvalue(fixed(8, 7; rounding=RoundNearest)(q)) == refraw(fixed(8, 7; rounding=RoundNearest), Rational{BigInt}(float(q)))
        @test FPN.Q0f15(fixed(16, 15)(0.25)) === FPN.Q0f15(0.25)
        @test FPN.Q0f7(fixed(16, 15; rounding=RoundNearest)(0.3)) === reinterpret(FPN.Q0f7, Int8(38))
        @test convert(fixed(16, 15), q) === fixed(16, 15)(q)
    end

    @testset "ambiguities" begin
        @test isempty(Test.detect_ambiguities(FixedPointArithmetic))
    end
end

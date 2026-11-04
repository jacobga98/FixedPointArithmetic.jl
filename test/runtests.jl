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

    @testset "no silent floating point" begin
        x = fixed(16, 15)(0.5)
        for f in (y -> y + 0.5, y -> 0.5 + y, y -> y * 0.5f0, y -> y < 0.1, y -> y == 0.5, y -> max(y, 0.1))
            @test_throws ArgumentError f(x)
        end
        for f in (sqrt, exp, log, sin, cos, inv, y -> y / y, y -> y / 2, y -> 2 / y, y -> y^3)
            @test_throws ArgumentError f(x)
        end
        z = Complex(x, x)
        @test_throws ArgumentError abs(z)
        @test_throws ArgumentError angle(z)
        @test_throws ArgumentError z / z
        @test_throws ArgumentError cis(x)
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

# FixedPointArithmetic.jl

Binary fixed-point numbers for simulating fixed-point signal processing in Julia: receivers, decoders, filters. A
format has a word length, a fraction length, a rounding mode and an overflow mode. Algorithms written generically
for floating point run unchanged on it, and arrays of a fixed-point type quantize whatever is stored into them.

**Status:** first draft, not registered. The API will change.

```julia
using FixedPointArithmetic

Q15  = fixed(16, 15)                                   # RoundDown, Saturate by default
LLR8 = fixed(8, 0; overflow=SaturateSymmetric)         # int8 soft bits in -127:127
LLR6 = fixed(6, 2; rounding=RoundNearest)              # 6-bit LLRs with step 0.25

a = Q15.(0.3 .* randn(1000))                           # quantize: round down, saturate
b = Q15.(0.3 .* randn(1000))
p = a .* b                                             # exact products, Fx{32, 30, …}
y = similar(a); y .= a .* b                            # stored back as Q15: rounded once
z = Complex(a[1], b[1]) * Complex(b[2], a[2])          # complex product: exact, 32-bit sums
Complex{Q15}(z)                                        # one rounding per component
```

## Arithmetic model

The model follows hardware fixed-point types such as `ap_fixed` and `ac_fixed`, and MATLAB's `fimath` with full
precision products and `KeepLSB` sums.

| Operation | Result |
|---|---|
| `X(x)`, `convert(X, x)`, storing into an array of `X` | `x` rounded with the rounding mode of `X`, then brought into its range with its overflow mode. Never throws: NaN becomes 0, and so do ±Inf with `Wrap`. |
| `a + b`, `a - b`, `-a`, `abs(a)` of one format | the same format, with its overflow mode. Loops like `acc += x` keep their type; choose the accumulator format by the type of `acc`. |
| `a * b` | exact, in a format with the word and fraction lengths added (at most 64 bits) |
| two formats with the same modes | promoted to the format that holds both exactly; the result keeps it |
| formats with different modes | no promotion: convert explicitly |
| `a + n`, `a * n`, `a < n` with an integer `n` | exact, in the format of `a` |
| `a + 0.5`, `a < 0.1`, … with a float | `ArgumentError`. Write the constant as `X(0.5)`, which also works for float `X`. |
| `/`, `sqrt`, `exp`, `log`, `sin`, `atan`, `abs(::Complex)`, `angle`, `cis`, `^p` | `ArgumentError`. Base would compute them in floating point; a fixed-point algorithm needs its own method. `x^2` works. |
| `ldexp(a, e)` | `a·2^e` in the format of `a`, with its rounding and overflow modes: the scaling step of block floating point |
| `float(a)`, `Float64(a)` | the exact value when the float has enough significand bits |

Rounding modes are the Base ones: `RoundDown` (floor, an arithmetic right shift), `RoundUp`, `RoundToZero`,
`RoundNearest` (ties to even), `RoundNearestTiesAway` and `RoundNearestTiesUp` (round half up). Overflow modes are
`Wrap`, `Saturate` and `SaturateSymmetric`.

The stored integer is available with `rawvalue(x)`, and `fromraw(X, i)` builds a value from one, for example to
compare with integer samples from another implementation.

## Performance and GPUs

Operations are thin wrappers around integer instructions and are branch-free on the data.

- **CPU vectorisation**, checked on Apple M1 with Julia 1.13:
  - a Q15 saturating add compiles to `sqadd`
  - a Q15 multiply stored back to Q15 compiles to `smull`, an arithmetic shift and `smin`
  - int16 → int8 symmetric saturation compiles to `smax`/`smin` and a narrow
- **GPU:** conversions never throw, so the type works inside KernelAbstractions kernels. Checked on Metal for
  complex products, saturating adds, float quantization and narrowing; all matched the CPU bit for bit. CUDA is not
  tested yet.

## Relation to FixedPointNumbers.jl

FixedPointNumbers' `Fixed{T,f}` is a storage format: its products keep the operand type, its conversions throw on
overflow, and its word length is the width of the storage type. This package converts to and from it when
FixedPointNumbers is loaded (`Fx(x::Fixed)`, `Fixed{T,f}(x::Fx)`), but the arithmetic is separate.

## Defaults chosen in this draft

These are open for discussion:
- **Default modes:** `fixed(W, F)` rounds down and saturates. Rounding down is what a right shift does; saturating
  is the safe overflow behaviour.
- **Sums keep their format** instead of growing by one bit as in `ap_fixed`, so that accumulation loops are
  type-stable.
- **Floats do not mix with fixed-point values.** Mixing would otherwise silently turn a fixed-point computation into a
  floating-point one.
- **Integers act as exact constants in the fixed-point format.** That keeps generic code such as `2x`, `-1 .* x`
  and `x > 0` working.

## Not yet

- Block-exponent helpers (an exponent per block from its maximum magnitude or mean power)
- Range and overflow logging, to choose formats from simulation data
- Non-power-of-two scaling (MATLAB's slope) and saturation limits below the format maximum, as in int8 LLRs
  limited to ±120
- Unsigned formats, and words longer than 64 bits
- Documentation build, CI, and GPU tests on CUDA

## Tests

```
julia --project -e 'using Pkg; Pkg.test()'
```

The tests compare every operation with an exact reference in `BigInt` and `Rational` arithmetic. They cover all
rounding and overflow modes, and all values of the 6- and 8-bit formats.

# FixedPointArithmetic.jl

Binary fixed-point numbers for simulating fixed-point signal processing in Julia: receivers, decoders, filters. A
format has a word length, a fraction length, a rounding mode and an overflow mode. Algorithms written generically
for floating point run unchanged on it, and arrays of a fixed-point type quantize whatever is stored into them.

!!! warning "Status"
    First draft, not registered. The API will change.

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

See [Arithmetic model](@ref) for what each operation returns, and the [API reference](@ref) for the functions.

## Performance and GPUs

Operations are thin wrappers around integer instructions and are branch-free on the data.

- **CPU vectorisation**, checked on Apple M5 with Julia 1.13:
  - a Q15 saturating add compiles to `sqadd`
  - a Q15 multiply stored back to Q15 compiles to `smull`, an arithmetic shift and `smin`
  - int16 → int8 symmetric saturation compiles to `smax`/`smin` and a narrow
- **GPU:** conversions never throw, so the type works inside KernelAbstractions kernels. Checked on Metal for
  complex products, saturating adds, float quantization, narrowing, division (including by zero), complex division
  and square roots; all matched the CPU bit for bit. CUDA is not tested yet.
- **Division and square roots** are integer algorithms: one integer division per quotient, and a bit-by-bit square
  root with one step per two bits of the operand. They are much slower than the other operations, as in hardware.

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
- **Quotients and square roots keep the format of the numerator or argument,** like sums, instead of growing as
  `ap_fixed` quotients do. [`divide`](@ref) takes another format explicitly.

## Not yet

- Block-exponent helpers (an exponent per block from its maximum magnitude or mean power)
- Range and overflow logging, to choose formats from simulation data
- Non-power-of-two scaling (MATLAB's slope) and saturation limits below the format maximum, as in int8 LLRs
  limited to ±120
- Unsigned formats, and words longer than 64 bits
- GPU tests on CUDA

## Tests

```
julia --project -e 'using Pkg; Pkg.test()'
```

The tests compare every operation with an exact reference in `BigInt` and `Rational` arithmetic. They cover all
rounding and overflow modes, and all values of the 6- and 8-bit formats.

# FixedPointArithmetic.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jacobga98.github.io/FixedPointArithmetic.jl/dev/)
[![CI](https://github.com/jacobga98/FixedPointArithmetic.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jacobga98/FixedPointArithmetic.jl/actions/workflows/CI.yml?query=branch%3Amain)

Binary fixed-point numbers for simulating fixed-point signal processing in Julia: receivers, decoders, filters. A
format has a word length, a fraction length, a rounding mode and an overflow mode. Algorithms written generically
for floating point run unchanged on it, and arrays of a fixed-point type quantize whatever is stored into them.

**Status:** first draft, not registered. The API will change.

```julia
using FixedPointArithmetic

Q15 = fixed(16, 15)          # RoundDown, Saturate by default
a = Q15.(0.3 .* randn(1000)) # quantize: round down, saturate
p = a .* a                   # exact products, Fx{32, 30, …}
y = similar(a); y .= p       # stored back as Q15: rounded once
```

See the [documentation](https://jacobga98.github.io/FixedPointArithmetic.jl/dev/) for the arithmetic model, performance
notes and the API reference.

# This file is a part of Julia. License is MIT: https://julialang.org/license

module Math

export sin, cos, sincos, tan, sinh, cosh, tanh, asin, acos, atan,
       asinh, acosh, atanh, sec, csc, cot, asec, acsc, acot,
       sech, csch, coth, asech, acsch, acoth,
       sinpi, cospi, sinc, cosc,
       cosd, cotd, cscd, secd, sind, tand,
       acosd, acotd, acscd, asecd, asind, atand,
       rad2deg, deg2rad,
       log, log2, log10, log1p, exponent, exp, exp2, exp10, expm1,
       cbrt, sqrt, significand,
       hypot, max, min, minmax, ldexp, frexp,
       clamp, clamp!, modf, ^, mod2pi, rem2pi,
       @evalpoly

import .Base: log, exp, sin, cos, tan, sinh, cosh, tanh, asin,
             acos, atan, asinh, acosh, atanh, sqrt, log2, log10,
             max, min, minmax, ^, exp2, muladd, rem,
             exp10, expm1, log1p

using .Base: sign_mask, exponent_mask, exponent_one,
            exponent_half, uinttype, significand_mask

using Core.Intrinsics: sqrt_llvm

using .Base: IEEEFloat

@noinline function throw_complex_domainerror(f::Symbol, x)
    throw(DomainError(x, string("$f will only return a complex result if called with a ",
                                "complex argument. Try $f(Complex(x)).")))
end
@noinline function throw_exp_domainerror(x)
    throw(DomainError(x, string("Exponentiation yielding a complex result requires a ",
                                "complex argument.\nReplace x^y with (x+0im)^y, ",
                                "Complex(x)^y, or similar.")))
end

for T in (Float16, Float32, Float64)
    @eval significand_bits(::Type{$T}) = $(trailing_ones(significand_mask(T)))
    @eval exponent_bits(::Type{$T}) = $(sizeof(T)*8 - significand_bits(T) - 1)
    @eval exponent_bias(::Type{$T}) = $(Int(exponent_one(T) >> significand_bits(T)))
    # maximum float exponent
    @eval exponent_max(::Type{$T}) = $(Int(exponent_mask(T) >> significand_bits(T)) - exponent_bias(T))
    # maximum float exponent without bias
    @eval exponent_raw_max(::Type{$T}) = $(Int(exponent_mask(T) >> significand_bits(T)))
end

# non-type specific math functions

"""
    clamp(x, lo, hi)

Return `x` if `lo <= x <= hi`. If `x > hi`, return `hi`. If `x < lo`, return `lo`. Arguments
are promoted to a common type.

# Examples
```jldoctest
julia> clamp.([pi, 1.0, big(10.)], 2., 9.)
3-element Array{BigFloat,1}:
 3.141592653589793238462643383279502884197169399375105820974944592307816406286198
 2.0
 9.0

julia> clamp.([11,8,5],10,6) # an example where lo > hi
3-element Array{Int64,1}:
  6
  6
 10
```
"""
clamp(x::X, lo::L, hi::H) where {X,L,H} =
    ifelse(x > hi, convert(promote_type(X,L,H), hi),
           ifelse(x < lo,
                  convert(promote_type(X,L,H), lo),
                  convert(promote_type(X,L,H), x)))

"""
    clamp!(array::AbstractArray, lo, hi)

Restrict values in `array` to the specified range, in-place.
See also [`clamp`](@ref).
"""
function clamp!(x::AbstractArray, lo, hi)
    @inbounds for i in eachindex(x)
        x[i] = clamp(x[i], lo, hi)
    end
    x
end

"""
    @horner(x, p...)
    Evaluate p[1] + x * (p[2] + x * (....)), i.e. a polynomial via Horner's rule
"""
macro horner(x, p...)
    ex = esc(p[end])
    for i = length(p)-1:-1:1
        ex = :(muladd(t, $ex, $(esc(p[i]))))
    end
    ex = quote local r = $ex end # structure this to add exactly one line number node for the macro
    return Expr(:block, :(local t = $(esc(x))), ex, :r)
end

# Evaluate p[1] + z*p[2] + z^2*p[3] + ... + z^(n-1)*p[n].  This uses
# Horner's method if z is real, but for complex z it uses a more
# efficient algorithm described in Knuth, TAOCP vol. 2, section 4.6.4,
# equation (3).

"""
    @evalpoly(z, c...)

Evaluate the polynomial ``\\sum_k c[k] z^{k-1}`` for the coefficients `c[1]`, `c[2]`, ...;
that is, the coefficients are given in ascending order by power of `z`.  This macro expands
to efficient inline code that uses either Horner's method or, for complex `z`, a more
efficient Goertzel-like algorithm.

# Examples
```jldoctest
julia> @evalpoly(3, 1, 0, 1)
10

julia> @evalpoly(2, 1, 0, 1)
5

julia> @evalpoly(2, 1, 1, 1)
7
```
"""
macro evalpoly(z, p...)
    a = :($(esc(p[end])))
    b = :($(esc(p[end-1])))
    as = []
    for i = length(p)-2:-1:1
        ai = Symbol("a", i)
        push!(as, :($ai = $a))
        a = :(muladd(r, $ai, $b))
        b = :($(esc(p[i])) - s * $ai) # see issue #15985 on fused mul-subtract
    end
    ai = :a0
    push!(as, :($ai = $a))
    C = Expr(:block,
             :(x = real(tt)),
             :(y = imag(tt)),
             :(r = x + x),
             :(s = muladd(x, x, y*y)),
             as...,
             :(muladd($ai, tt, $b)))
    R = Expr(:macrocall, Symbol("@horner"), (), :tt, map(esc, p)...)
    :(let tt = $(esc(z))
          isa(tt, Complex) ? $C : $R
      end)
end

"""
    rad2deg(x)

Convert `x` from radians to degrees.

# Examples
```jldoctest
julia> rad2deg(pi)
180.0
```
"""
rad2deg(z::AbstractFloat) = z * (180 / oftype(z, pi))

"""
    deg2rad(x)

Convert `x` from degrees to radians.

# Examples
```jldoctest
julia> deg2rad(90)
1.5707963267948966
```
"""
deg2rad(z::AbstractFloat) = z * (oftype(z, pi) / 180)
rad2deg(z::Real) = rad2deg(float(z))
deg2rad(z::Real) = deg2rad(float(z))
rad2deg(z::Number) = (z/pi)*180
deg2rad(z::Number) = (z*pi)/180

log(b::T, x::T) where {T<:Number} = log(x)/log(b)

"""
    log(b,x)

Compute the base `b` logarithm of `x`. Throws [`DomainError`](@ref) for negative
[`Real`](@ref) arguments.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> log(4,8)
1.5

julia> log(4,2)
0.5

julia> log(-2, 3)
ERROR: DomainError with -2.0:
log will only return a complex result if called with a complex argument. Try log(Complex(x)).
Stacktrace:
 [1] throw_complex_domainerror(::Symbol, ::Float64) at ./math.jl:31
[...]

julia> log(2, -3)
ERROR: DomainError with -3.0:
log will only return a complex result if called with a complex argument. Try log(Complex(x)).
Stacktrace:
 [1] throw_complex_domainerror(::Symbol, ::Float64) at ./math.jl:31
[...]
```

!!! note
    If `b` is a power of 2 or 10, [`log2`](@ref) or [`log10`](@ref) should be used, as these will
    typically be faster and more accurate. For example,

    ```jldoctest
    julia> log(100,1000000)
    2.9999999999999996

    julia> log10(1000000)/2
    3.0
    ```
"""
log(b::Number, x::Number) = log(promote(b,x)...)

# type specific math functions

const libm = Base.libm_name

# functions with no domain error
"""
    sinh(x)

Compute hyperbolic sine of `x`.
"""
sinh(x::Number)

"""
    cosh(x)

Compute hyperbolic cosine of `x`.
"""
cosh(x::Number)

"""
    tanh(x)

Compute hyperbolic tangent of `x`.
"""
tanh(x::Number)

"""
    atan(y)
    atan(y, x)

Compute the inverse tangent of `y` or `y/x`, respectively.

For one argument, this is the angle in radians between the positive *x*-axis and the point
(1, *y*), returning a value in the interval ``[-\\pi/2, \\pi/2]``.

For two arguments, this is the angle in radians between the positive *x*-axis and the
point (*x*, *y*), returning a value in the interval ``[-\\pi, \\pi]``. This corresponds to a
standard [`atan2`](https://en.wikipedia.org/wiki/Atan2) function.
"""
atan(x::Number)

"""
    asinh(x)

Compute the inverse hyperbolic sine of `x`.
"""
asinh(x::Number)

"""
    expm1(x)

Accurately compute ``e^x-1``.
"""
expm1(x)
for f in (:exp2, :expm1)
    @eval begin
        ($f)(x::Float64) = ccall(($(string(f)),libm), Float64, (Float64,), x)
        ($f)(x::Float32) = ccall(($(string(f,"f")),libm), Float32, (Float32,), x)
        ($f)(x::Real) = ($f)(float(x))
    end
end

"""
    exp2(x)

Compute the base 2 exponential of `x`, in other words ``2^x``.

# Examples
```jldoctest
julia> exp2(5)
32.0
```
"""
exp2(x::AbstractFloat) = 2^x

"""
    exp10(x)

Compute the base 10 exponential of `x`, in other words ``10^x``.

# Examples
```jldoctest
julia> exp10(2)
100.0
```
"""
exp10(x::AbstractFloat) = 10^x

for f in (:sinh, :cosh, :tanh, :atan, :asinh, :exp, :expm1)
    @eval ($f)(x::AbstractFloat) = error("not implemented for ", typeof(x))
end

# functions with special cases for integer arguments
@inline function exp2(x::Base.BitInteger)
    if x > 1023
        Inf64
    elseif x <= -1023
        # if -1073 < x <= -1023 then Result will be a subnormal number
        # Hex literal with padding must be used to work on 32bit machine
        reinterpret(Float64, 0x0000_0000_0000_0001 << ((x + 1074) % UInt))
    else
        # We will cast everything to Int64 to avoid errors in case of Int128
        # If x is a Int128, and is outside the range of Int64, then it is not -1023<x<=1023
        reinterpret(Float64, (exponent_bias(Float64) + (x % Int64)) << (significand_bits(Float64) % UInt))
    end
end

# utility for converting NaN return to DomainError
# the branch in nan_dom_err prevents its callers from inlining, so be sure to force it
# until the heuristics can be improved
@inline nan_dom_err(out, x) = isnan(out) & !isnan(x) ? throw(DomainError(x, "NaN result for non-NaN input.")) : out

# functions that return NaN on non-NaN argument for domain error
"""
    sin(x)

Compute sine of `x`, where `x` is in radians.
"""
sin(x::Number)

"""
    cos(x)

Compute cosine of `x`, where `x` is in radians.
"""
cos(x::Number)

"""
    tan(x)

Compute tangent of `x`, where `x` is in radians.
"""
tan(x::Number)

"""
    asin(x)

Compute the inverse sine of `x`, where the output is in radians.
"""
asin(x::Number)

"""
    acos(x)

Compute the inverse cosine of `x`, where the output is in radians
"""
acos(x::Number)

"""
    acosh(x)

Compute the inverse hyperbolic cosine of `x`.
"""
acosh(x::Number)

"""
    atanh(x)

Compute the inverse hyperbolic tangent of `x`.
"""
atanh(x::Number)

"""
    log(x)

Compute the natural logarithm of `x`. Throws [`DomainError`](@ref) for negative
[`Real`](@ref) arguments. Use complex negative arguments to obtain complex results.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> log(2)
0.6931471805599453

julia> log(-3)
ERROR: DomainError with -3.0:
log will only return a complex result if called with a complex argument. Try log(Complex(x)).
Stacktrace:
 [1] throw_complex_domainerror(::Symbol, ::Float64) at ./math.jl:31
[...]
```
"""
log(x::Number)

"""
    log2(x)

Compute the logarithm of `x` to base 2. Throws [`DomainError`](@ref) for negative
[`Real`](@ref) arguments.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> log2(4)
2.0

julia> log2(10)
3.321928094887362

julia> log2(-2)
ERROR: DomainError with -2.0:
NaN result for non-NaN input.
Stacktrace:
 [1] nan_dom_err at ./math.jl:325 [inlined]
[...]
```
"""
log2(x)

"""
    log10(x)

Compute the logarithm of `x` to base 10.
Throws [`DomainError`](@ref) for negative [`Real`](@ref) arguments.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> log10(100)
2.0

julia> log10(2)
0.3010299956639812

julia> log10(-2)
ERROR: DomainError with -2.0:
NaN result for non-NaN input.
Stacktrace:
 [1] nan_dom_err at ./math.jl:325 [inlined]
[...]
```
"""
log10(x)

"""
    log1p(x)

Accurate natural logarithm of `1+x`. Throws [`DomainError`](@ref) for [`Real`](@ref)
arguments less than -1.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> log1p(-0.5)
-0.6931471805599453

julia> log1p(0)
0.0

julia> log1p(-2)
ERROR: DomainError with -2.0:
log1p will only return a complex result if called with a complex argument. Try log1p(Complex(x)).
Stacktrace:
 [1] throw_complex_domainerror(::Symbol, ::Float64) at ./math.jl:31
[...]
```
"""
log1p(x)
for f in (:log2, :log10)
    @eval begin
        @inline ($f)(x::Float64) = nan_dom_err(ccall(($(string(f)), libm), Float64, (Float64,), x), x)
        @inline ($f)(x::Float32) = nan_dom_err(ccall(($(string(f, "f")), libm), Float32, (Float32,), x), x)
        @inline ($f)(x::Real) = ($f)(float(x))
    end
end

@inline function sqrt(x::Union{Float32,Float64})
    x < zero(x) && throw_complex_domainerror(:sqrt, x)
    sqrt_llvm(x)
end

"""
    sqrt(x)

Return ``\\sqrt{x}``. Throws [`DomainError`](@ref) for negative [`Real`](@ref) arguments.
Use complex negative arguments instead. The prefix operator `√` is equivalent to `sqrt`.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> sqrt(big(81))
9.0

julia> sqrt(big(-81))
ERROR: DomainError with -81.0:
NaN result for non-NaN input.
Stacktrace:
 [1] sqrt(::BigFloat) at ./mpfr.jl:501
[...]

julia> sqrt(big(complex(-81)))
0.0 + 9.0im
```
"""
sqrt(x::Real) = sqrt(float(x))

"""
    hypot(x, y)

Compute the hypotenuse ``\\sqrt{x^2+y^2}`` avoiding overflow and underflow.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> a = 10^10;

julia> hypot(a, a)
1.4142135623730951e10

julia> √(a^2 + a^2) # a^2 overflows
ERROR: DomainError with -2.914184810805068e18:
sqrt will only return a complex result if called with a complex argument. Try sqrt(Complex(x)).
Stacktrace:
[...]
```
"""
hypot(x::Number, y::Number) = hypot(promote(x, y)...)
function hypot(x::T, y::T) where T<:Number
    ax = abs(x)
    ay = abs(y)
    if ax < ay
        ax, ay = ay, ax
    end
    if iszero(ax)
        r = ay / oneunit(ax)
    else
        r = ay / ax
    end

    rr = ax * sqrt(1 + r * r)

    # Use type of rr to make sure that return type is the same for
    # all branches
    if isnan(r)
        isinf(ax) && return oftype(rr, Inf)
        isinf(ay) && return oftype(rr, Inf)
        return oftype(rr, r)
    else
        return rr
    end
end

"""
    hypot(x...)

Compute the hypotenuse ``\\sqrt{\\sum x_i^2}`` avoiding overflow and underflow.
"""
hypot(x::Number...) = sqrt(sum(abs2(y) for y in x))

atan(y::Real, x::Real) = atan(promote(float(y),float(x))...)
atan(y::T, x::T) where {T<:AbstractFloat} = Base.no_op_err("atan", T)

max(x::T, y::T) where {T<:AbstractFloat} = ifelse((y > x) | (signbit(y) < signbit(x)),
                                    ifelse(isnan(x), x, y), ifelse(isnan(y), y, x))


min(x::T, y::T) where {T<:AbstractFloat} = ifelse((y < x) | (signbit(y) > signbit(x)),
                                    ifelse(isnan(x), x, y), ifelse(isnan(y), y, x))

minmax(x::T, y::T) where {T<:AbstractFloat} =
    ifelse(isnan(x) | isnan(y), ifelse(isnan(x), (x,x), (y,y)),
           ifelse((y > x) | (signbit(x) > signbit(y)), (x,y), (y,x)))


"""
    ldexp(x, n)

Compute ``x \\times 2^n``.

# Examples
```jldoctest
julia> ldexp(5., 2)
20.0
```
"""
function ldexp(x::T, e::Integer) where T<:IEEEFloat
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if k == 0 # x is subnormal
        xs == 0 && return x # +-0
        m = leading_zeros(xs) - exponent_bits(T)
        ys = xs << unsigned(m)
        xu = ys | (xu & sign_mask(T))
        k = 1 - m
        # underflow, otherwise may have integer underflow in the following n + k
        e < -50000 && return flipsign(T(0.0), x)
    end
    # For cases where e of an Integer larger than Int make sure we properly
    # overflow/underflow; this is optimized away otherwise.
    if e > typemax(Int)
        return flipsign(T(Inf), x)
    elseif e < typemin(Int)
        return flipsign(T(0.0), x)
    end
    n = e % Int
    k += n
    # overflow, if k is larger than maximum possible exponent
    if k >= exponent_raw_max(T)
        return flipsign(T(Inf), x)
    end
    if k > 0 # normal case
        xu = (xu & ~exponent_mask(T)) | (rem(k, uinttype(T)) << significand_bits(T))
        return reinterpret(T, xu)
    else # subnormal case
        if k <= -significand_bits(T) # underflow
            # overflow, for the case of integer overflow in n + k
            e > 50000 && return flipsign(T(Inf), x)
            return flipsign(T(0.0), x)
        end
        k += significand_bits(T)
        z = T(2.0)^-significand_bits(T)
        xu = (xu & ~exponent_mask(T)) | (rem(k, uinttype(T)) << significand_bits(T))
        return z*reinterpret(T, xu)
    end
end
ldexp(x::Float16, q::Integer) = Float16(ldexp(Float32(x), q))

"""
    exponent(x) -> Int

Get the exponent of a normalized floating-point number.
"""
function exponent(x::T) where T<:IEEEFloat
    @noinline throw1(x) = throw(DomainError(x, "Cannot be NaN or Inf."))
    @noinline throw2(x) = throw(DomainError(x, "Cannot be subnormal converted to 0."))
    xs = reinterpret(Unsigned, x) & ~sign_mask(T)
    xs >= exponent_mask(T) && throw1(x)
    k = Int(xs >> significand_bits(T))
    if k == 0 # x is subnormal
        xs == 0 && throw2(x)
        m = leading_zeros(xs) - exponent_bits(T)
        k = 1 - m
    end
    return k - exponent_bias(T)
end

"""
    significand(x)

Extract the `significand(s)` (a.k.a. mantissa), in binary representation, of a
floating-point number. If `x` is a non-zero finite number, then the result will be
a number of the same type on the interval ``[1,2)``. Otherwise `x` is returned.

# Examples
```jldoctest
julia> significand(15.2)/15.2
0.125

julia> significand(15.2)*8
15.2
```
"""
function significand(x::T) where T<:IEEEFloat
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x # NaN or Inf
    if xs <= (~exponent_mask(T) & ~sign_mask(T)) # x is subnormal
        xs == 0 && return x # +-0
        m = unsigned(leading_zeros(xs) - exponent_bits(T))
        xs <<= m
        xu = xs | (xu & sign_mask(T))
    end
    xu = (xu & ~exponent_mask(T)) | exponent_one(T)
    return reinterpret(T, xu)
end

"""
    frexp(val)

Return `(x,exp)` such that `x` has a magnitude in the interval ``[1/2, 1)`` or 0,
and `val` is equal to ``x \\times 2^{exp}``.
"""
function frexp(x::T) where T<:IEEEFloat
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x, 0 # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if k == 0 # x is subnormal
        xs == 0 && return x, 0 # +-0
        m = leading_zeros(xs) - exponent_bits(T)
        xs <<= unsigned(m)
        xu = xs | (xu & sign_mask(T))
        k = 1 - m
    end
    k -= (exponent_bias(T) - 1)
    xu = (xu & ~exponent_mask(T)) | exponent_half(T)
    return reinterpret(T, xu), k
end

"""
    rem(x, y, r::RoundingMode)

Compute the remainder of `x` after integer division by `y`, with the quotient rounded
according to the rounding mode `r`. In other words, the quantity

    x - y*round(x/y,r)

without any intermediate rounding.

- if `r == RoundNearest`, then the result is exact, and in the interval
  ``[-|y|/2, |y|/2]``. See also [`RoundNearest`](@ref).

- if `r == RoundToZero` (default), then the result is exact, and in the interval
  ``[0, |y|)`` if `x` is positive, or ``(-|y|, 0]`` otherwise. See also [`RoundToZero`](@ref).

- if `r == RoundDown`, then the result is in the interval ``[0, y)`` if `y` is positive, or
  ``(y, 0]`` otherwise. The result may not be exact if `x` and `y` have different signs, and
  `abs(x) < abs(y)`. See also[`RoundDown`](@ref).

- if `r == RoundUp`, then the result is in the interval `(-y,0]` if `y` is positive, or
  `[0,-y)` otherwise. The result may not be exact if `x` and `y` have the same sign, and
  `abs(x) < abs(y)`. See also [`RoundUp`](@ref).

"""
rem(x, y, ::RoundingMode{:ToZero}) = rem(x,y)
rem(x, y, ::RoundingMode{:Down}) = mod(x,y)
rem(x, y, ::RoundingMode{:Up}) = mod(x,-y)

rem(x::Float64, y::Float64, ::RoundingMode{:Nearest}) =
    ccall((:remainder, libm),Float64,(Float64,Float64),x,y)
rem(x::Float32, y::Float32, ::RoundingMode{:Nearest}) =
    ccall((:remainderf, libm),Float32,(Float32,Float32),x,y)
rem(x::Float16, y::Float16, r::RoundingMode{:Nearest}) = Float16(rem(Float32(x), Float32(y), r))


"""
    modf(x)

Return a tuple `(fpart, ipart)` of the fractional and integral parts of a number. Both parts
have the same sign as the argument.

# Examples
```jldoctest
julia> modf(3.5)
(0.5, 3.0)

julia> modf(-3.5)
(-0.5, -3.0)
```
"""
modf(x) = rem(x,one(x)), trunc(x)

function modf(x::Float32)
    temp = Ref{Float32}()
    f = ccall((:modff, libm), Float32, (Float32, Ptr{Float32}), x, temp)
    f, temp[]
end

function modf(x::Float64)
    temp = Ref{Float64}()
    f = ccall((:modf, libm), Float64, (Float64, Ptr{Float64}), x, temp)
    f, temp[]
end

@inline function ^(x::Float64, y::Float64)
    z = ccall("llvm.pow.f64", llvmcall, Float64, (Float64, Float64), x, y)
    if isnan(z) & !isnan(x+y)
        throw_exp_domainerror(x)
    end
    z
end
@inline function ^(x::Float32, y::Float32)
    z = ccall("llvm.pow.f32", llvmcall, Float32, (Float32, Float32), x, y)
    if isnan(z) & !isnan(x+y)
        throw_exp_domainerror(x)
    end
    z
end
@inline ^(x::Float64, y::Integer) = ccall("llvm.pow.f64", llvmcall, Float64, (Float64, Float64), x, Float64(y))
@inline ^(x::Float32, y::Integer) = ccall("llvm.pow.f32", llvmcall, Float32, (Float32, Float32), x, Float32(y))
@inline ^(x::Float16, y::Integer) = Float16(Float32(x) ^ y)
@inline literal_pow(::typeof(^), x::Float16, ::Val{p}) where {p} = Float16(literal_pow(^,Float32(x),Val(p)))

## rem2pi-related calculations ##

function add22condh(xh::Float64, xl::Float64, yh::Float64, yl::Float64)
    # This algorithm, due to Dekker, computes the sum of two
    # double-double numbers and returns the high double. References:
    # [1] http://www.digizeitschriften.de/en/dms/img/?PID=GDZPPN001170007
    # [2] https://doi.org/10.1007/BF01397083
    r = xh+yh
    s = (abs(xh) > abs(yh)) ? (xh-r+yh+yl+xl) : (yh-r+xh+xl+yl)
    zh = r+s
    return zh
end

# multiples of pi/2, as double-double (ie with "tail")
const pi1o2_h  = 1.5707963267948966     # convert(Float64, pi * BigFloat(1/2))
const pi1o2_l  = 6.123233995736766e-17  # convert(Float64, pi * BigFloat(1/2) - pi1o2_h)

const pi2o2_h  = 3.141592653589793      # convert(Float64, pi * BigFloat(1))
const pi2o2_l  = 1.2246467991473532e-16 # convert(Float64, pi * BigFloat(1) - pi2o2_h)

const pi3o2_h  = 4.71238898038469       # convert(Float64, pi * BigFloat(3/2))
const pi3o2_l  = 1.8369701987210297e-16 # convert(Float64, pi * BigFloat(3/2) - pi3o2_h)

const pi4o2_h  = 6.283185307179586      # convert(Float64, pi * BigFloat(2))
const pi4o2_l  = 2.4492935982947064e-16 # convert(Float64, pi * BigFloat(2) - pi4o2_h)

"""
    rem2pi(x, r::RoundingMode)

Compute the remainder of `x` after integer division by `2π`, with the quotient rounded
according to the rounding mode `r`. In other words, the quantity

    x - 2π*round(x/(2π),r)

without any intermediate rounding. This internally uses a high precision approximation of
2π, and so will give a more accurate result than `rem(x,2π,r)`

- if `r == RoundNearest`, then the result is in the interval ``[-π, π]``. This will generally
  be the most accurate result. See also [`RoundNearest`](@ref).

- if `r == RoundToZero`, then the result is in the interval ``[0, 2π]`` if `x` is positive,.
  or ``[-2π, 0]`` otherwise. See also [`RoundToZero`](@ref).

- if `r == RoundDown`, then the result is in the interval ``[0, 2π]``.
  See also [`RoundDown`](@ref).
- if `r == RoundUp`, then the result is in the interval ``[-2π, 0]``.
  See also [`RoundUp`](@ref).

# Examples
```jldoctest
julia> rem2pi(7pi/4, RoundNearest)
-0.7853981633974485

julia> rem2pi(7pi/4, RoundDown)
5.497787143782138
```
"""
function rem2pi end
function rem2pi(x::Float64, ::RoundingMode{:Nearest})
    abs(x) < pi && return x

    n,y = rem_pio2_kernel(x)

    if iseven(n)
        if n & 2 == 2 # n % 4 == 2: add/subtract pi
            if y.hi <= 0
                return add22condh(y.hi,y.lo,pi2o2_h,pi2o2_l)
            else
                return add22condh(y.hi,y.lo,-pi2o2_h,-pi2o2_l)
            end
        else          # n % 4 == 0: add 0
            return y.hi+y.lo
        end
    else
        if n & 2 == 2 # n % 4 == 3: subtract pi/2
            return add22condh(y.hi,y.lo,-pi1o2_h,-pi1o2_l)
        else          # n % 4 == 1: add pi/2
            return add22condh(y.hi,y.lo,pi1o2_h,pi1o2_l)
        end
    end
end
function rem2pi(x::Float64, ::RoundingMode{:ToZero})
    ax = abs(x)
    ax <= 2*Float64(pi,RoundDown) && return x

    n,y = rem_pio2_kernel(x)

    if iseven(n)
        if n & 2 == 2 # n % 4 == 2: add pi
            z = add22condh(y.hi,y.lo,pi2o2_h,pi2o2_l)
        else          # n % 4 == 0: add 0 or 2pi
            if y.hi > 0
                z = y.hi+y.lo
            else      # negative: add 2pi
                z = add22condh(y.hi,y.lo,pi4o2_h,pi4o2_l)
            end
        end
    else
        if n & 2 == 2 # n % 4 == 3: add 3pi/2
            z = add22condh(y.hi,y.lo,pi3o2_h,pi3o2_l)
        else          # n % 4 == 1: add pi/2
            z = add22condh(y.hi,y.lo,pi1o2_h,pi1o2_l)
        end
    end
    copysign(z,x)
end
function rem2pi(x::Float64, ::RoundingMode{:Down})
    if x < pi4o2_h
        if x >= 0
            return x
        elseif x > -pi4o2_h
            return add22condh(x,0.0,pi4o2_h,pi4o2_l)
        end
    end

    n,y = rem_pio2_kernel(x)

    if iseven(n)
        if n & 2 == 2 # n % 4 == 2: add pi
            return add22condh(y.hi,y.lo,pi2o2_h,pi2o2_l)
        else          # n % 4 == 0: add 0 or 2pi
            if y.hi > 0
                return y.hi+y.lo
            else      # negative: add 2pi
                return add22condh(y.hi,y.lo,pi4o2_h,pi4o2_l)
            end
        end
    else
        if n & 2 == 2 # n % 4 == 3: add 3pi/2
            return add22condh(y.hi,y.lo,pi3o2_h,pi3o2_l)
        else          # n % 4 == 1: add pi/2
            return add22condh(y.hi,y.lo,pi1o2_h,pi1o2_l)
        end
    end
end
function rem2pi(x::Float64, ::RoundingMode{:Up})
    if x > -pi4o2_h
        if x <= 0
            return x
        elseif x < pi4o2_h
            return add22condh(x,0.0,-pi4o2_h,-pi4o2_l)
        end
    end

    n,y = rem_pio2_kernel(x)

    if iseven(n)
        if n & 2 == 2 # n % 4 == 2: sub pi
            return add22condh(y.hi,y.lo,-pi2o2_h,-pi2o2_l)
        else          # n % 4 == 0: sub 0 or 2pi
            if y.hi < 0
                return y.hi+y.lo
            else      # positive: sub 2pi
                return add22condh(y.hi,y.lo,-pi4o2_h,-pi4o2_l)
            end
        end
    else
        if n & 2 == 2 # n % 4 == 3: sub pi/2
            return add22condh(y.hi,y.lo,-pi1o2_h,-pi1o2_l)
        else          # n % 4 == 1: sub 3pi/2
            return add22condh(y.hi,y.lo,-pi3o2_h,-pi3o2_l)
        end
    end
end

rem2pi(x::Float32, r::RoundingMode) = Float32(rem2pi(Float64(x), r))
rem2pi(x::Float16, r::RoundingMode) = Float16(rem2pi(Float64(x), r))
rem2pi(x::Int32, r::RoundingMode) = rem2pi(Float64(x), r)
function rem2pi(x::Int64, r::RoundingMode)
    fx = Float64(x)
    fx == x || throw(ArgumentError("Int64 argument to rem2pi is too large: $x"))
    rem2pi(fx, r)
end

"""
    mod2pi(x)

Modulus after division by `2π`, returning in the range ``[0,2π)``.

This function computes a floating point representation of the modulus after division by
numerically exact `2π`, and is therefore not exactly the same as `mod(x,2π)`, which would
compute the modulus of `x` relative to division by the floating-point number `2π`.

# Examples
```jldoctest
julia> mod2pi(9*pi/4)
0.7853981633974481
```
"""
mod2pi(x) = rem2pi(x,RoundDown)

# generic fallback; for number types, promotion.jl does promotion

"""
    muladd(x, y, z)

Combined multiply-add: computes `x*y+z`, but allowing the add and multiply to be merged
with each other or with surrounding operations for performance.
For example, this may be implemented as an [`fma`](@ref) if the hardware supports it
efficiently.
The result can be different on different machines and can also be different on the same machine
due to constant propagation or other optimizations.
See [`fma`](@ref).

# Examples
```jldoctest
julia> muladd(3, 2, 1)
7

julia> 3 * 2 + 1
7
```
"""
muladd(x,y,z) = x*y+z

# Float16 definitions

for func in (:sin,:cos,:tan,:asin,:acos,:atan,:sinh,:cosh,:tanh,:asinh,:acosh,
             :atanh,:exp,:exp2,:exp10,:log,:log2,:log10,:sqrt,:lgamma,:log1p)
    @eval begin
        $func(a::Float16) = Float16($func(Float32(a)))
        $func(a::ComplexF16) = ComplexF16($func(ComplexF32(a)))
    end
end

for func in (:atan,:hypot)
    @eval begin
        $func(a::Float16,b::Float16) = Float16($func(Float32(a),Float32(b)))
    end
end

cbrt(a::Float16) = Float16(cbrt(Float32(a)))
sincos(a::Float16) = Float16.(sincos(Float32(a)))

# helper functions for Libm functionality

"""
    highword(x)

Return the high word of `x` as a `UInt32`.
"""
@inline highword(x::Float64) = highword(reinterpret(UInt64, x))
@inline highword(x::UInt64)  = (x >>> 32) % UInt32
@inline highword(x::Float32) = reinterpret(UInt32, x)

@inline fromhighword(::Type{Float64}, u::UInt32) = reinterpret(Float64, UInt64(u) << 32)
@inline fromhighword(::Type{Float32}, u::UInt32) = reinterpret(Float32, u)


"""
    poshighword(x)

Return positive part of the high word of `x` as a `UInt32`.
"""
@inline poshighword(x::Float64) = poshighword(reinterpret(UInt64, x))
@inline poshighword(x::UInt64)  = highword(x) & 0x7fffffff
@inline poshighword(x::Float32) = highword(x) & 0x7fffffff

# More special functions
include("special/cbrt.jl")
include("special/exp.jl")
include("special/exp10.jl")
include("special/ldexp_exp.jl")
include("special/hyperbolic.jl")
include("special/trig.jl")
include("special/rem_pio2.jl")
include("special/log.jl")

# `missing` definitions for functions in this module
for f in (:(acos), :(acosh), :(asin), :(asinh), :(atan), :(atanh),
          :(sin), :(sinh), :(cos), :(cosh), :(tan), :(tanh),
          :(exp), :(exp2), :(expm1), :(log), :(log10), :(log1p),
          :(log2), :(exponent), :(sqrt))
    @eval $(f)(::Missing) = missing
end

end # module


# This file is a part of Julia. License is MIT: https://julialang.org/license


# This type must be kept in sync with the C struct in src/gc.h
struct GC_Num
    allocd          ::Int64 # GC internal
    deferred_alloc  ::Int64 # GC internal
    freed           ::Int64 # GC internal
    malloc          ::Int64
    realloc         ::Int64
    poolalloc       ::Int64
    bigalloc        ::Int64
    freecall        ::Int64
    total_time      ::Int64
    total_allocd    ::Int64 # GC internal
    since_sweep     ::Int64 # GC internal
    collect         ::Csize_t # GC internal
    pause           ::Cint
    full_sweep      ::Cint
end

gc_num() = ccall(:jl_gc_num, GC_Num, ())

# This type is to represent differences in the counters, so fields may be negative
struct GC_Diff
    allocd      ::Int64 # Bytes allocated
    malloc      ::Int64 # Number of GC aware malloc()
    realloc     ::Int64 # Number of GC aware realloc()
    poolalloc   ::Int64 # Number of pool allocation
    bigalloc    ::Int64 # Number of big (non-pool) allocation
    freecall    ::Int64 # Number of GC aware free()
    total_time  ::Int64 # Time spent in garbage collection
    pause       ::Int64 # Number of GC pauses
    full_sweep  ::Int64 # Number of GC full collection
end

gc_total_bytes(gc_num::GC_Num) =
    gc_num.allocd + gc_num.deferred_alloc + gc_num.total_allocd

function GC_Diff(new::GC_Num, old::GC_Num)
    # logic from `src/gc.c:jl_gc_total_bytes`
    old_allocd = gc_total_bytes(old)
    new_allocd = gc_total_bytes(new)
    return GC_Diff(new_allocd - old_allocd,
                   new.malloc       - old.malloc,
                   new.realloc      - old.realloc,
                   new.poolalloc    - old.poolalloc,
                   new.bigalloc     - old.bigalloc,
                   new.freecall     - old.freecall,
                   new.total_time   - old.total_time,
                   new.pause        - old.pause,
                   new.full_sweep   - old.full_sweep)
end

function gc_alloc_count(diff::GC_Diff)
    diff.malloc + diff.realloc + diff.poolalloc + diff.bigalloc
end


# total time spend in garbage collection, in nanoseconds
gc_time_ns() = ccall(:jl_gc_total_hrtime, UInt64, ())

"""
    Base.gc_live_bytes()

Return the total size (in bytes) of objects currently in memory.
This is computed as the total size of live objects after
the last garbage collection, plus the number of bytes allocated
since then.
"""
function gc_live_bytes()
    num = gc_num()
    Int(ccall(:jl_gc_live_bytes, Int64, ())) + num.allocd + num.deferred_alloc
end

# print elapsed time, return expression value
const _mem_units = ["byte", "KiB", "MiB", "GiB", "TiB", "PiB"]
const _cnt_units = ["", " k", " M", " G", " T", " P"]
function prettyprint_getunits(value, numunits, factor)
    if value == 0 || value == 1
        return (value, 1)
    end
    unit = ceil(Int, log(value) / log(factor))
    unit = min(numunits, unit)
    number = value/factor^(unit-1)
    return number, unit
end

function padded_nonzero_print(value, str)
    if value != 0
        blanks = "                "[1:(18 - length(str))]
        println(str, ":", blanks, value)
    end
end


function format_bytes(bytes) # also used by InteractiveUtils
    bytes, mb = prettyprint_getunits(bytes, length(_mem_units), Int64(1024))
    if mb == 1
        return string(Int(bytes), " ", _mem_units[mb], bytes==1 ? "" : "s")
    else
        return string(Ryu.writefixed(Float64(bytes), 3), " ", _mem_units[mb])
    end
end

function time_print(elapsedtime, bytes=0, gctime=0, allocs=0)
    timestr = Ryu.writefixed(Float64(elapsedtime/1e9), 6)
    length(timestr) < 10 && print(" "^(10 - length(timestr)))
    print(timestr, " seconds")
    if bytes != 0 || allocs != 0
        allocs, ma = prettyprint_getunits(allocs, length(_cnt_units), Int64(1000))
        if ma == 1
            print(" (", Int(allocs), _cnt_units[ma], allocs==1 ? " allocation: " : " allocations: ")
        else
            print(" (", Ryu.writefixed(Float64(allocs), 2), _cnt_units[ma], " allocations: ")
        end
        print(format_bytes(bytes))
    end
    if gctime > 0
        print(", ", Ryu.writefixed(Float64(100*gctime/elapsedtime), 2), "% gc time")
    end
    if bytes != 0 || allocs != 0
        print(")")
    end
    nothing
end

function timev_print(elapsedtime, diff::GC_Diff)
    allocs = gc_alloc_count(diff)
    time_print(elapsedtime, diff.allocd, diff.total_time, allocs)
    print("\nelapsed time (ns): $elapsedtime\n")
    padded_nonzero_print(diff.total_time,   "gc time (ns)")
    padded_nonzero_print(diff.allocd,       "bytes allocated")
    padded_nonzero_print(diff.poolalloc,    "pool allocs")
    padded_nonzero_print(diff.bigalloc,     "non-pool GC allocs")
    padded_nonzero_print(diff.malloc,       "malloc() calls")
    padded_nonzero_print(diff.realloc,      "realloc() calls")
    padded_nonzero_print(diff.freecall,     "free() calls")
    padded_nonzero_print(diff.pause,        "GC pauses")
    padded_nonzero_print(diff.full_sweep,   "full collections")
end

"""
    @time

A macro to execute an expression, printing the time it took to execute, the number of
allocations, and the total number of bytes its execution caused to be allocated, before
returning the value of the expression.

See also [`@timev`](@ref), [`@timed`](@ref), [`@elapsed`](@ref), and
[`@allocated`](@ref).

!!! note
    For more serious benchmarking, consider the `@btime` macro from the BenchmarkTools.jl
    package which among other things evaluates the function multiple times in order to
    reduce noise.

```julia-repl
julia> @time rand(10^6);
  0.001525 seconds (7 allocations: 7.630 MiB)

julia> @time begin
           sleep(0.3)
           1+1
       end
  0.301395 seconds (8 allocations: 336 bytes)
2
```
"""
macro time(ex)
    quote
        while false; end # compiler heuristic: compile this block (alter this if the heuristic changes)
        local stats = gc_num()
        local elapsedtime = time_ns()
        local val = $(esc(ex))
        elapsedtime = time_ns() - elapsedtime
        local diff = GC_Diff(gc_num(), stats)
        time_print(elapsedtime, diff.allocd, diff.total_time,
                   gc_alloc_count(diff))
        println()
        val
    end
end

"""
    @timev

This is a verbose version of the `@time` macro. It first prints the same information as
`@time`, then any non-zero memory allocation counters, and then returns the value of the
expression.

See also [`@time`](@ref), [`@timed`](@ref), [`@elapsed`](@ref), and
[`@allocated`](@ref).

```julia-repl
julia> @timev rand(10^6);
  0.001006 seconds (7 allocations: 7.630 MiB)
elapsed time (ns): 1005567
bytes allocated:   8000256
pool allocs:       6
malloc() calls:    1
```
"""
macro timev(ex)
    quote
        while false; end # compiler heuristic: compile this block (alter this if the heuristic changes)
        local stats = gc_num()
        local elapsedtime = time_ns()
        local val = $(esc(ex))
        elapsedtime = time_ns() - elapsedtime
        timev_print(elapsedtime, GC_Diff(gc_num(), stats))
        val
    end
end

"""
    @elapsed

A macro to evaluate an expression, discarding the resulting value, instead returning the
number of seconds it took to execute as a floating-point number.

See also [`@time`](@ref), [`@timev`](@ref), [`@timed`](@ref),
and [`@allocated`](@ref).

```julia-repl
julia> @elapsed sleep(0.3)
0.301391426
```
"""
macro elapsed(ex)
    quote
        while false; end # compiler heuristic: compile this block (alter this if the heuristic changes)
        local t0 = time_ns()
        $(esc(ex))
        (time_ns() - t0) / 1e9
    end
end

# total number of bytes allocated so far
gc_bytes(b::Ref{Int64}) = ccall(:jl_gc_get_total_bytes, Cvoid, (Ptr{Int64},), b)
# NOTE: gc_bytes() is deprecated
function gc_bytes()
    b = Ref{Int64}()
    gc_bytes(b)
    b[]
end

"""
    @allocated

A macro to evaluate an expression, discarding the resulting value, instead returning the
total number of bytes allocated during evaluation of the expression.

See also [`@time`](@ref), [`@timev`](@ref), [`@timed`](@ref),
and [`@elapsed`](@ref).

```julia-repl
julia> @allocated rand(10^6)
8000080
```
"""
macro allocated(ex)
    quote
        while false; end # compiler heuristic: compile this block (alter this if the heuristic changes)
        local b0 = Ref{Int64}(0)
        local b1 = Ref{Int64}(0)
        gc_bytes(b0)
        $(esc(ex))
        gc_bytes(b1)
        b1[] - b0[]
    end
end

"""
    @timed

A macro to execute an expression, and return the value of the expression, elapsed time,
total bytes allocated, garbage collection time, and an object with various memory allocation
counters.

See also [`@time`](@ref), [`@timev`](@ref), [`@elapsed`](@ref), and
[`@allocated`](@ref).

```julia-repl
julia> val, t, bytes, gctime, memallocs = @timed rand(10^6);

julia> t
0.006634834

julia> bytes
8000256

julia> gctime
0.0055765

julia> fieldnames(typeof(memallocs))
(:allocd, :malloc, :realloc, :poolalloc, :bigalloc, :freecall, :total_time, :pause, :full_sweep)

julia> memallocs.total_time
5576500
```
"""
macro timed(ex)
    quote
        while false; end # compiler heuristic: compile this block (alter this if the heuristic changes)
        local stats = gc_num()
        local elapsedtime = time_ns()
        local val = $(esc(ex))
        elapsedtime = time_ns() - elapsedtime
        local diff = GC_Diff(gc_num(), stats)
        val, elapsedtime/1e9, diff.allocd, diff.total_time/1e9, diff
    end
end


## printing with color ##

const text_colors = Dict{Union{Symbol,Int},String}(
    :black         => "\033[30m",
    :red           => "\033[31m",
    :green         => "\033[32m",
    :yellow        => "\033[33m",
    :blue          => "\033[34m",
    :magenta       => "\033[35m",
    :cyan          => "\033[36m",
    :white         => "\033[37m",
    :light_black   => "\033[90m", # gray
    :light_red     => "\033[91m",
    :light_green   => "\033[92m",
    :light_yellow  => "\033[93m",
    :light_blue    => "\033[94m",
    :light_magenta => "\033[95m",
    :light_cyan    => "\033[96m",
    :normal        => "\033[0m",
    :default       => "\033[39m",
    :bold          => "\033[1m",
    :underline     => "\033[4m",
    :blink         => "\033[5m",
    :reverse       => "\033[7m",
    :hidden        => "\033[8m",
    :nothing       => "",
)

for i in 0:255
    text_colors[i] = "\033[38;5;$(i)m"
end

const disable_text_style = Dict{Symbol,String}(
    :bold      => "\033[22m",
    :underline => "\033[24m",
    :blink     => "\033[25m",
    :reverse   => "\033[27m",
    :hidden    => "\033[28m",
    :normal    => "",
    :default   => "",
    :nothing   => "",
)

# Create a docstring with an automatically generated list
# of colors.
available_text_colors = collect(Iterators.filter(x -> !isa(x, Integer), keys(text_colors)))
const possible_formatting_symbols = [:normal, :bold, :default]
available_text_colors = cat(
    sort!(intersect(available_text_colors, possible_formatting_symbols), rev=true),
    sort!(setdiff(  available_text_colors, possible_formatting_symbols));
    dims=1)

const available_text_colors_docstring =
    string(join([string("`:", key,"`")
                 for key in available_text_colors], ",\n", ", or \n"))

"""Dictionary of color codes for the terminal.

Available colors are: $available_text_colors_docstring as well as the integers 0 to 255 inclusive.

The color `:default` will print text in the default color while the color `:normal`
will print text with all text properties (like boldness) reset.
Printing with the color `:nothing` will print the string without modifications.
"""
text_colors

function with_output_color(f::Function, color::Union{Int, Symbol}, io::IO, args...; bold::Bool = false)
    buf = IOBuffer()
    iscolor = get(io, :color, false)
    try f(IOContext(buf, io), args...)
    finally
        str = String(take!(buf))
        if !iscolor
            print(io, str)
        else
            bold && color === :bold && (color = :nothing)
            enable_ansi  = get(text_colors, color, text_colors[:default]) *
                               (bold ? text_colors[:bold] : "")
            disable_ansi = (bold ? disable_text_style[:bold] : "") *
                               get(disable_text_style, color, text_colors[:default])
            first = true
            for line in split(str, '\n')
                first || print(buf, '\n')
                first = false
                isempty(line) && continue
                print(buf, enable_ansi, line, disable_ansi)
            end
            print(io, String(take!(buf)))
        end
    end
end

"""
    printstyled([io], xs...; bold::Bool=false, color::Union{Symbol,Int}=:normal)

Print `xs` in a color specified as a symbol or integer, optionally in bold.

`color` may take any of the values $(Base.available_text_colors_docstring)
or an integer between 0 and 255 inclusive. Note that not all terminals support 256 colors.
If the keyword `bold` is given as `true`, the result will be printed in bold.
"""
printstyled(io::IO, msg...; bold::Bool=false, color::Union{Int,Symbol}=:normal) =
    with_output_color(print, color, io, msg...; bold=bold)
printstyled(msg...; bold::Bool=false, color::Union{Int,Symbol}=:normal) =
    printstyled(stdout, msg...; bold=bold, color=color)

"""
    Base.julia_cmd(juliapath=joinpath(Sys.BINDIR::String, julia_exename()))

Return a julia command similar to the one of the running process.
Propagates any of the `--cpu-target`, `--sysimage`, `--compile`, `--sysimage-native-code`,
`--compiled-modules`, `--inline`, `--check-bounds`, `--optimize`, `-g`,
`--code-coverage`, and `--depwarn`
command line arguments that are not at their default values.

Among others, `--math-mode`, `--warn-overwrite`, and `--trace-compile` are notably not propagated currently.

!!! compat "Julia 1.1"
    Only the `--cpu-target`, `--sysimage`, `--depwarn`, `--compile` and `--check-bounds` flags were propagated before Julia 1.1.
"""
function julia_cmd(julia=joinpath(Sys.BINDIR::String, julia_exename()))
    opts = JLOptions()
    cpu_target = unsafe_string(opts.cpu_target)
    image_file = unsafe_string(opts.image_file)
    addflags = String[]
    let compile = if opts.compile_enabled == 0
                      "no"
                  elseif opts.compile_enabled == 2
                      "all"
                  elseif opts.compile_enabled == 3
                      "min"
                  else
                      "" # default = "yes"
                  end
        isempty(compile) || push!(addflags, "--compile=$compile")
    end
    let depwarn = if opts.depwarn == 0
                      "no"
                  elseif opts.depwarn == 2
                      "error"
                  else
                      "" # default = "yes"
                  end
        isempty(depwarn) || push!(addflags, "--depwarn=$depwarn")
    end
    let check_bounds = if opts.check_bounds == 1
                      "yes" # on
                  elseif opts.check_bounds == 2
                      "no" # off
                  else
                      "" # "default"
                  end
        isempty(check_bounds) || push!(addflags, "--check-bounds=$check_bounds")
    end
    opts.can_inline == 0 && push!(addflags, "--inline=no")
    opts.use_compiled_modules == 0 && push!(addflags, "--compiled-modules=no")
    opts.opt_level == 2 || push!(addflags, "-O$(opts.opt_level)")
    push!(addflags, "-g$(opts.debug_level)")
    if opts.code_coverage != 0
        # Forward the code-coverage flag only if applicable (if the filename is pid-dependent)
        coverage_file = (opts.output_code_coverage != C_NULL) ?  unsafe_string(opts.output_code_coverage) : ""
        if isempty(coverage_file) || occursin("%p", coverage_file)
            if opts.code_coverage == 1
                push!(addflags, "--code-coverage=user")
            elseif opts.code_coverage == 2
                push!(addflags, "--code-coverage=all")
            end
            isempty(coverage_file) || push!(addflags, "--code-coverage=$coverage_file")
        end
    end
    if opts.malloc_log != 0
        if opts.malloc_log == 1
            push!(addflags, "--track-allocation=user")
        elseif opts.malloc_log == 2
            push!(addflags, "--track-allocation=all")
        end
    end
    return `$julia -C$cpu_target -J$image_file $addflags`
end

function julia_exename()
    if ccall(:jl_is_debugbuild, Cint, ()) == 0
        return @static Sys.iswindows() ? "julia.exe" : "julia"
    else
        return @static Sys.iswindows() ? "julia-debug.exe" : "julia-debug"
    end
end

"""
    securezero!(o)

`securezero!` fills the memory associated with an object `o` with zeros.
Unlike `fill!(o,0)` and similar code, which might be optimized away by
the compiler for objects about to be discarded, the `securezero!` function
will always be called.
"""
function securezero! end
@noinline securezero!(a::AbstractArray{<:Number}) = fill!(a, 0)
@noinline unsafe_securezero!(p::Ptr{T}, len::Integer=1) where {T} =
    ccall(:memset, Ptr{T}, (Ptr{T}, Cint, Csize_t), p, 0, len*sizeof(T))
unsafe_securezero!(p::Ptr{Cvoid}, len::Integer=1) = Ptr{Cvoid}(unsafe_securezero!(Ptr{UInt8}(p), len))

"""
    Base.getpass(message::AbstractString) -> Base.SecretBuffer

Display a message and wait for the user to input a secret, returning an `IO`
object containing the secret.

Note that on Windows, the secret might be displayed as it is typed; see
`Base.winprompt` for securely retrieving username/password pairs from a
graphical interface.
"""
function getpass end

if Sys.iswindows()
function getpass(input::TTY, output::IO, prompt::AbstractString)
    input === stdin || throw(ArgumentError("getpass only works for stdin"))
    print(output, prompt, ": ")
    flush(output)
    s = SecretBuffer()
    plen = 0
    while true
        c = UInt8(ccall(:_getch, Cint, ()))
        if c == 0xff || c == UInt8('\n') || c == UInt8('\r')
            break # EOF or return
        elseif c == 0x00 || c == 0xe0
            ccall(:_getch, Cint, ()) # ignore function/arrow keys
        elseif c == UInt8('\b') && plen > 0
            plen -= 1 # delete last character on backspace
        elseif !iscntrl(Char(c)) && plen < 128
            write(s, c)
        end
    end
    return seekstart(s)
end
else
function getpass(input::TTY, output::IO, prompt::AbstractString)
    (input === stdin && output === stdout) || throw(ArgumentError("getpass only works for stdin"))
    msg = string(prompt, ": ")
    unsafe_SecretBuffer!(ccall(:getpass, Cstring, (Cstring,), msg))
end
end

# allow new getpass methods to be defined if stdin has been
# redirected to some custom stream, e.g. in IJulia.
getpass(prompt::AbstractString) = getpass(stdin, stdout, prompt)

"""
    prompt(message; default="") -> Union{String, Nothing}

Displays the `message` then waits for user input. Input is terminated when a newline (\\n)
is encountered or EOF (^D) character is entered on a blank line. If a `default` is provided
then the user can enter just a newline character to select the `default`.

See also `Base.getpass` and `Base.winprompt` for secure entry of passwords.
"""
function prompt(input::IO, output::IO, message::AbstractString; default::AbstractString="")
    msg = !isempty(default) ? "$message [$default]: " : "$message: "
    print(output, msg)
    uinput = readline(input, keep=true)
    isempty(uinput) && return nothing  # Encountered an EOF
    uinput = chomp(uinput)
    isempty(uinput) ? default : uinput
end

# allow new prompt methods to be defined if stdin has been
# redirected to some custom stream, e.g. in IJulia.
prompt(message::AbstractString; default::AbstractString="") = prompt(stdin, stdout, message, default=default)

# Windows authentication prompt
if Sys.iswindows()
    struct CREDUI_INFO
        cbSize::UInt32
        parent::Ptr{Cvoid}
        pszMessageText::Ptr{UInt16}
        pszCaptionText::Ptr{UInt16}
        banner::Ptr{Cvoid}
    end

    const CREDUIWIN_GENERIC                 = 0x0001
    const CREDUIWIN_IN_CRED_ONLY            = 0x0020
    const CREDUIWIN_ENUMERATE_CURRENT_USER  = 0x0200

    const CRED_PACK_GENERIC_CREDENTIALS     = 0x0004

    const ERROR_SUCCESS                     = 0x0000
    const ERROR_CANCELLED                   = 0x04c7

    function winprompt(message, caption, default_username; prompt_username = true)
        # Step 1: Create an encrypted username/password bundle that will be used to set
        #         the default username (in theory could also provide a default password)
        credbuf = Vector{UInt8}(undef, 1024)
        credbufsize = Ref{UInt32}(sizeof(credbuf))
        succeeded = ccall((:CredPackAuthenticationBufferW, "credui.dll"), stdcall, Bool,
            (UInt32, Cwstring, Cwstring, Ptr{UInt8}, Ptr{UInt32}),
             CRED_PACK_GENERIC_CREDENTIALS, default_username, "", credbuf, credbufsize)
        @assert succeeded

        # Step 2: Create the actual dialog
        #      2.1: Set up the window
        messageArr = Base.cwstring(message)
        captionArr = Base.cwstring(caption)
        pfSave = Ref{Bool}(false)
        cred = Ref{CREDUI_INFO}(CREDUI_INFO(sizeof(CREDUI_INFO), C_NULL, pointer(messageArr), pointer(captionArr), C_NULL))
        dwflags = CREDUIWIN_GENERIC | CREDUIWIN_ENUMERATE_CURRENT_USER
        if !prompt_username
            # Disable setting anything other than default_username
            dwflags |= CREDUIWIN_IN_CRED_ONLY
        end
        authPackage = Ref{Culong}(0)
        outbuf_data = Ref{Ptr{Cvoid}}(C_NULL)
        outbuf_size = Ref{Culong}(0)

        #      2.2: Do the actual request
        code = ccall((:CredUIPromptForWindowsCredentialsW, "credui.dll"), stdcall, UInt32, (Ptr{CREDUI_INFO}, UInt32, Ptr{Culong},
            Ptr{Cvoid}, Culong, Ptr{Ptr{Cvoid}}, Ptr{Culong}, Ptr{Bool}, UInt32), cred, 0, authPackage, credbuf, credbufsize[],
            outbuf_data, outbuf_size, pfSave, dwflags)

        #      2.3: If that failed for any reason other than the user canceling, error out.
        #           If the user canceled, just return nothing
        code == ERROR_CANCELLED && return nothing
        windowserror(:winprompt, code != ERROR_SUCCESS)

        # Step 3: Convert encrypted credentials back to plain text
        passbuf = Vector{UInt16}(undef, 1024)
        passlen = Ref{UInt32}(length(passbuf))
        usernamebuf = Vector{UInt16}(undef, 1024)
        usernamelen = Ref{UInt32}(length(usernamebuf))
        # Need valid buffers for domain, even though we don't care
        dummybuf = Vector{UInt16}(undef, 1024)
        succeeded = ccall((:CredUnPackAuthenticationBufferW, "credui.dll"), Bool,
            (UInt32, Ptr{Cvoid}, UInt32, Ptr{UInt16}, Ptr{UInt32}, Ptr{UInt16}, Ptr{UInt32}, Ptr{UInt16}, Ptr{UInt32}),
            0, outbuf_data[], outbuf_size[], usernamebuf, usernamelen, dummybuf, Ref{UInt32}(1024), passbuf, passlen)
        windowserror(:winprompt, !succeeded)

        # Step 4: Free the encrypted buffer
        # ccall(:SecureZeroMemory, Ptr{Cvoid}, (Ptr{Cvoid}, Csize_t), outbuf_data[], outbuf_size[]) - not an actual function
        unsafe_securezero!(outbuf_data[], outbuf_size[])
        ccall((:CoTaskMemFree, "ole32.dll"), Cvoid, (Ptr{Cvoid},), outbuf_data[])

        # Done.
        passbuf_ = passbuf[1:passlen[]-1]
        result = (String(transcode(UInt8, usernamebuf[1:usernamelen[]-1])),
                  SecretBuffer!(transcode(UInt8, passbuf_)))
        securezero!(passbuf_)
        securezero!(passbuf)

        return result
    end

end

unsafe_crc32c(a, n, crc) = ccall(:jl_crc32c, UInt32, (UInt32, Ptr{UInt8}, Csize_t), crc, a, n)

_crc32c(a::Union{Array{UInt8},FastContiguousSubArray{UInt8,N,<:Array{UInt8}} where N}, crc::UInt32=0x00000000) =
    unsafe_crc32c(a, length(a) % Csize_t, crc)

_crc32c(s::String, crc::UInt32=0x00000000) = unsafe_crc32c(s, sizeof(s) % Csize_t, crc)

function _crc32c(io::IO, nb::Integer, crc::UInt32=0x00000000)
    nb < 0 && throw(ArgumentError("number of bytes to checksum must be ≥ 0, got $nb"))
    # use block size 24576=8192*3, since that is the threshold for
    # 3-way parallel SIMD code in the underlying jl_crc32c C function.
    buf = Vector{UInt8}(undef, min(nb, 24576))
    while !eof(io) && nb > 24576
        n = readbytes!(io, buf)
        crc = unsafe_crc32c(buf, n, crc)
        nb -= n
    end
    return unsafe_crc32c(buf, readbytes!(io, buf, min(nb, length(buf))), crc)
end
_crc32c(io::IO, crc::UInt32=0x00000000) = _crc32c(io, typemax(Int64), crc)
_crc32c(io::IOStream, crc::UInt32=0x00000000) = _crc32c(io, filesize(io)-position(io), crc)
_crc32c(uuid::UUID, crc::UInt32=0x00000000) =
    ccall(:jl_crc32c, UInt32, (UInt32, Ref{UInt128}, Csize_t), crc, uuid.value, 16)

"""
    @kwdef typedef

This is a helper macro that automatically defines a keyword-based constructor for the type
declared in the expression `typedef`, which must be a `struct` or `mutable struct`
expression. The default argument is supplied by declaring fields of the form `field::T =
default` or `field = default`. If no default is provided then the keyword argument becomes
a required keyword argument in the resulting type constructor.

Inner constructors can still be defined, but at least one should accept arguments in the
same form as the default inner constructor (i.e. one positional argument per field) in
order to function correctly with the keyword outer constructor.

!!! compat "Julia 1.1"
    `Base.@kwdef` for parametric structs, and structs with supertypes
    requires at least Julia 1.1.

# Examples
```jldoctest
julia> Base.@kwdef struct Foo
           a::Int = 1         # specified default
           b::String          # required keyword
       end
Foo

julia> Foo(b="hi")
Foo(1, "hi")

julia> Foo()
ERROR: UndefKeywordError: keyword argument b not assigned
Stacktrace:
[...]
```
"""
macro kwdef(expr)
    expr = macroexpand(__module__, expr) # to expand @static
    expr isa Expr && expr.head === :struct || error("Invalid usage of @kwdef")
    T = expr.args[2]
    if T isa Expr && T.head === :<:
        T = T.args[1]
    end

    params_ex = Expr(:parameters)
    call_args = Any[]

    _kwdef!(expr.args[3], params_ex.args, call_args)
    # Only define a constructor if the type has fields, otherwise we'll get a stack
    # overflow on construction
    if !isempty(params_ex.args)
        if T isa Symbol
            kwdefs = :(($(esc(T)))($params_ex) = ($(esc(T)))($(call_args...)))
        elseif T isa Expr && T.head === :curly
            # if T == S{A<:AA,B<:BB}, define two methods
            #   S(...) = ...
            #   S{A,B}(...) where {A<:AA,B<:BB} = ...
            S = T.args[1]
            P = T.args[2:end]
            Q = [U isa Expr && U.head === :<: ? U.args[1] : U for U in P]
            SQ = :($S{$(Q...)})
            kwdefs = quote
                ($(esc(S)))($params_ex) =($(esc(S)))($(call_args...))
                ($(esc(SQ)))($params_ex) where {$(esc.(P)...)} =
                    ($(esc(SQ)))($(call_args...))
            end
        else
            error("Invalid usage of @kwdef")
        end
    else
        kwdefs = nothing
    end
    quote
        Base.@__doc__($(esc(expr)))
        $kwdefs
    end
end

# @kwdef helper function
# mutates arguments inplace
function _kwdef!(blk, params_args, call_args)
    for i in eachindex(blk.args)
        ei = blk.args[i]
        if ei isa Symbol
            #  var
            push!(params_args, ei)
            push!(call_args, ei)
        elseif ei isa Expr
            if ei.head === :(=)
                lhs = ei.args[1]
                if lhs isa Symbol
                    #  var = defexpr
                    var = lhs
                elseif lhs isa Expr && lhs.head === :(::) && lhs.args[1] isa Symbol
                    #  var::T = defexpr
                    var = lhs.args[1]
                else
                    # something else, e.g. inline inner constructor
                    #   F(...) = ...
                    continue
                end
                defexpr = ei.args[2]  # defexpr
                push!(params_args, Expr(:kw, var, esc(defexpr)))
                push!(call_args, var)
                blk.args[i] = lhs
            elseif ei.head === :(::) && ei.args[1] isa Symbol
                # var::Typ
                var = ei.args[1]
                push!(params_args, var)
                push!(call_args, var)
            elseif ei.head === :block
                # can arise with use of @static inside type decl
                _kwdef!(ei, params_args, call_args)
            end
        end
    end
    blk
end

# testing

"""
    Base.runtests(tests=["all"]; ncores=ceil(Int, Sys.CPU_THREADS / 2),
                  exit_on_error=false, [seed])

Run the Julia unit tests listed in `tests`, which can be either a string or an array of
strings, using `ncores` processors. If `exit_on_error` is `false`, when one test
fails, all remaining tests in other files will still be run; they are otherwise discarded,
when `exit_on_error == true`.
If a seed is provided via the keyword argument, it is used to seed the
global RNG in the context where the tests are run; otherwise the seed is chosen randomly.
"""
function runtests(tests = ["all"]; ncores = ceil(Int, Sys.CPU_THREADS / 2),
                  exit_on_error=false,
                  seed::Union{BitInteger,Nothing}=nothing)
    if isa(tests,AbstractString)
        tests = split(tests)
    end
    exit_on_error && push!(tests, "--exit-on-error")
    seed !== nothing && push!(tests, "--seed=0x$(string(seed % UInt128, base=16))") # cast to UInt128 to avoid a minus sign
    ENV2 = copy(ENV)
    ENV2["JULIA_CPU_THREADS"] = "$ncores"
    try
        run(setenv(`$(julia_cmd()) $(joinpath(Sys.BINDIR::String,
            Base.DATAROOTDIR, "julia", "test", "runtests.jl")) $tests`, ENV2))
        nothing
    catch
        buf = PipeBuffer()
        Base.require(Base, :InteractiveUtils).versioninfo(buf)
        error("A test has failed. Please submit a bug report (https://github.com/JuliaLang/julia/issues)\n" *
              "including error messages above and the output of versioninfo():\n$(read(buf, String))")
    end
end


##############
# Reflection #
##############

"""
    Cassette.Reflection

A struct representing the information retrieved via `Cassette.reflect`.

A `Reflection` is essentially just a convenient bundle of information about a
specific method invocation.

# Fields

- `signature`: the invocation signature (in `Tuple{...}` type form) for the invoked method.

- `method`: the `Method` object associated with the invoked method.

- `static_params`: a `Vector` representing the invoked method's static parameter list.

- `code_info`: the `CodeInfo` object associated with the invoked method.
"""
mutable struct Reflection
    signature::DataType
    method::Method
    static_params::Vector{Any}
    code_info::CodeInfo
end

if VERSION < v"1.1.0-DEV.762"
    copy_code_info(code_info) = Core.Compiler.copy_code_info(code_info)
else
    copy_code_info(code_info) = copy(code_info)
end

if VERSION < v"1.2.0-DEV.573"
    specialize_method(method, metharg, methsp, world, force) = Core.Compiler.code_for_method(method, metharg, methsp, world, force)
else
    specialize_method(method, metharg, methsp, world, force) = Core.Compiler.specialize_method(method, metharg, methsp, force)
end



# Return `Reflection` for signature `sigtypes` and `world`, if possible. Otherwise, return `nothing`.
function reflect(@nospecialize(sigtypes::Tuple), world::UInt = typemax(UInt))
    if length(sigtypes) > 2 && sigtypes[1] === typeof(invoke)
        @assert sigtypes[3] <: Type{<:Tuple}
        sigtypes = (sigtypes[2], sigtypes[3].parameters[1].parameters...)
    end
    # This works around a subtyping bug. Basically, callers can deconstruct upstream
    # `UnionAll` types in such a way that results in a type with free type variables, in
    # which case subtyping can just break.
    #
    # God help you if you try to use a type parameter here (e.g. `::Type{S} where S<:Tuple`)
    # instead of this nutty workaround, because the compiler can just rewrite `S` into
    # whatever it thinks is "type equal" to the actual provided value. In other words, if
    # `S` is defined as e.g. `f(::Type{S}) where S`, and you call `f(T)`, you should NOT
    # assume that `S === T`. If you did, SHAME ON YOU. It doesn't matter that such an
    # assumption holds true for essentially all other kinds of values. I haven't counted in
    # a while, but I'm pretty sure I have ~40+ hellish years of Julia experience, and this
    # still catches me every time. Who even uses this crazy language?
    S = Tuple{map(s -> Core.Compiler.has_free_typevars(s) ? typeof(s.parameters[1]) : s, sigtypes)...}
    (S.parameters[1]::DataType).name.module === Core.Compiler && return nothing
    _methods = Base._methods_by_ftype(S, -1, world)
    method_index = 0
    for i in 1:length(_methods)
        if _methods[i][1] === S
            method_index = i
            break
        end
    end
    method_index === 0 && return nothing
    type_signature, raw_static_params, method = _methods[method_index]
    method_instance = specialize_method(method, type_signature, raw_static_params, world, false)
    method_instance === nothing && return nothing
    method_signature = method.sig
    static_params = Any[raw_static_params...]
    code_info = Core.Compiler.retrieve_code_info(method_instance)
    isa(code_info, CodeInfo) || return nothing
    code_info = copy_code_info(code_info)
@static if VERSION >= v"1.3.0-DEV.379"
        edges = Core.MethodInstance[method_instance]
        @assert code_info.edges === nothing
        code_info.edges = edges
    end
    return Reflection(S, method, static_params, code_info)
end

###########
# overdub #
###########

"""
    Cassette.OVERDUB_CONTEXT_NAME

The variable name bound to `overdub`'s `Context` argument in its `@generated`
method definition.

This binding can be used to manually reference/destructure `overdub` arguments
within `Expr` thunks emitted by user-provided passes.

See also: [`OVERDUB_ARGUMENTS_NAME`](@ref), [`@pass`](@ref), [`overdub`](@ref)
"""
const OVERDUB_CONTEXT_NAME = gensym("overdub_context")

"""
    Cassette.OVERDUB_ARGUMENTS_NAME

The variable name bound to `overdub`'s tuple of non-`Context` arguments in its
`@generated` method definition.

This binding can be used to manually reference/destructure `overdub` arguments
within `Expr` thunks emitted by user-provided passes.

See also: [`OVERDUB_CONTEXT_NAME`](@ref), [`@pass`](@ref), [`overdub`](@ref)
"""
const OVERDUB_ARGUMENTS_NAME = gensym("overdub_arguments")

# The `overdub` pass has four intertwined tasks:
#   1. Apply the user-provided pass, if one is given
#   2. Munge the reflection-generated IR into a valid form for returning from
#      `overdub_generator` (i.e. add new argument slots, substitute static
#      parameters, destructure overdub arguments into underlying method slots, etc.)
#   3. Perform the statement replacement central to the overdubbing pass (see `overdub` docstring)
#   4. If tagging is enabled, do the necessary IR transforms for the metadata tagging system
function overdub_pass!(reflection::Reflection,
                       context_type::DataType,
                       is_invoke::Bool = false)
    signature = reflection.signature
    method = reflection.method
    static_params = reflection.static_params
    code_info = reflection.code_info

    # TODO: This `iskwfunc` is part of a hack that `overdub_pass!` implements in order to fix
    # jrevels/Cassette.jl#48. The assumptions made by this hack are quite fragile, so we
    # should eventually get Base to expose a standard/documented API for this. Here, we see
    # this hack's first assumption: that `Core.kwfunc(f)` is going to return a function whose
    # type name is prefixed by `#kw##`. More assumptions for this hack will be commented on
    # as we go.
    name = String(signature.parameters[1].name.name)
    iskwfunc = startswith(name, "#kw##")
    istaggingenabled = hastagging(context_type)

    is_calloverload = false
    if length(code_info.slotnames) > 0
        # TODO: is this sufficient?
        is_calloverload = first(code_info.slotnames) !== Symbol("#self#")
    end

    #=== execute user-provided pass (is a no-op by default) ===#

    if !iskwfunc
        code_info = passtype(context_type)(context_type, reflection)
        isa(code_info, Expr) && return code_info
    end

    #=== munge the code into a valid form for `overdub_generator` ===#

    # NOTE: The slotflags set by this pass are set according to what makes sense based on the
    # compiler's actual `@code_lowered` output in practice, since this real-world output does
    # not seem to match Julia's developer documentation.

    # construct new slotnames/slotflags for added slots
    code_info.slotnames = Any[:overdub, OVERDUB_CONTEXT_NAME, OVERDUB_ARGUMENTS_NAME, code_info.slotnames...]
    code_info.slotflags = UInt8[0x00, 0x00, 0x00, code_info.slotflags...]
    n_prepended_slots = 3
    overdub_ctx_slot = SlotNumber(2)
    overdub_args_slot = SlotNumber(3)
    invoke_offset = is_invoke ? 2 : 0

    # For the sake of convenience, the rest of this pass will translate `code_info`'s fields
    # into these overdubbed equivalents instead of updating `code_info` in-place. Then, at
    # the end of the pass, we'll reset `code_info` fields accordingly.
    overdubbed_code = Any[]
    overdubbed_codelocs = Int32[]

    # destructure the generated argument slots into the overdubbed method's argument slots.
    n_actual_args = fieldcount(signature)
    n_method_args = Int(method.nargs)
    for i in 1:n_method_args
        slot = i + n_prepended_slots
        offset = i + invoke_offset
        if is_invoke && is_calloverload && i == 1
            offset = 2 # 1 is invoke, 2 is f, 3 is Tuple{}, 4... is args
        end
        actual_argument = Expr(:call, GlobalRef(Core, :getfield), overdub_args_slot, offset)
        push!(overdubbed_code, :($(SlotNumber(slot)) = $actual_argument))
        push!(overdubbed_codelocs, code_info.codelocs[1])
        code_info.slotflags[slot] |= 0x02 # ensure this slotflag has the "assigned" bit set
    end

    # If `method` is a varargs method, we have to restructure the original method call's
    # trailing arguments into a tuple and assign that tuple to the expected argument slot.
    if method.isva
        if !isempty(overdubbed_code)
            # remove the final slot reassignment leftover from the previous destructuring
            pop!(overdubbed_code)
            pop!(overdubbed_codelocs)
        end
        if hastagging(context_type)
            trailing_arguments = Expr(:call, GlobalRef(Cassette, :_tagged_new_tuple_unsafe), overdub_ctx_slot)
        else
            trailing_arguments = Expr(:call, GlobalRef(Core, :tuple))
        end
        for i in n_method_args:n_actual_args
            push!(overdubbed_code, Expr(:call, GlobalRef(Core, :getfield), overdub_args_slot, i + invoke_offset))
            push!(overdubbed_codelocs, code_info.codelocs[1])
            push!(trailing_arguments.args, SSAValue(length(overdubbed_code)))
        end
        push!(overdubbed_code, Expr(:(=), SlotNumber(n_method_args + n_prepended_slots), trailing_arguments))
        push!(overdubbed_codelocs, code_info.codelocs[1])
    end

    #=== finish initialization of `overdubbed_code`/`overdubbed_codelocs` ===#

    # substitute static parameters, offset slot numbers by number of added slots, and
    # offset statement indices by the number of additional statements
    Base.Meta.partially_inline!(code_info.code, Any[], method.sig, static_params,
                                n_prepended_slots, length(overdubbed_code), :propagate)

    original_code_start_index = length(overdubbed_code) + 1

    append!(overdubbed_code, code_info.code)
    append!(overdubbed_codelocs, code_info.codelocs)

    #=== perform tagged module transformation if tagging is enabled ===#

    if istaggingenabled && !iskwfunc
        # find `GlobalRef`s in IR and set up tagged module replacements
        modules = Any[]
        original_code_region = view(overdubbed_code, original_code_start_index:length(overdubbed_code))
        replace_match!(x -> isa(x, GlobalRef), original_code_region) do x
            x.mod === Core && in(x.name, (:tuple, :_apply)) && return x
            m = GlobalRef(parentmodule(x.mod), nameof(x.mod))
            i = findfirst(isequal(m), modules)
            if isa(i, Nothing)
                push!(modules, m)
                i = length(modules)
            end
            return Expr(:replaceglobalref, original_code_start_index + i - 1, x)
        end
        for i in 1:length(modules)
            modules[i] = Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, :fetch_tagged_module)), overdub_ctx_slot, modules[i])
        end

        # insert `fetch_tagged_module`s at the `original_code_start_index`
        insert_statements!(overdubbed_code, overdubbed_codelocs,
                            (x, i) -> i == original_code_start_index ? length(modules) + 1 : nothing,
                            (x, i) -> [modules..., x])

        # append `tagged_globalref_set_meta!` to `GlobalRef() = ...` statements
        insert_statements!(overdubbed_code, overdubbed_codelocs,
                            (x, i) -> begin
                                if (i > original_code_start_index &&
                                    Base.Meta.isexpr(x, :(=)) &&
                                    Base.Meta.isexpr(x.args[1], :replaceglobalref))
                                    return 3
                                end
                                return nothing
                            end,
                            (x, i) -> begin
                                lhs, rhs = x.args
                                tagmodssa = SSAValue(lhs.args[1])
                                globalref = lhs.args[2]
                                name = QuoteNode(globalref.name)
                                return [
                                    rhs,
                                    Expr(:(=), globalref, Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, :untag)), SSAValue(i), overdub_ctx_slot)),
                                    Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, :tagged_globalref_set_meta!)), overdub_ctx_slot, tagmodssa, name, SSAValue(i))
                                ]
                            end)

        # replace `GlobalRef`-loads with `Cassette.tagged_globalref`
        original_code_start_index += length(modules)
        insert_statements!(overdubbed_code, overdubbed_codelocs,
                            (x, i) -> begin
                                i >= original_code_start_index || return nothing
                                stmt = Base.Meta.isexpr(x, :(=)) ? x.args[2] : x
                                Base.Meta.isexpr(stmt, :replaceglobalref) && return 1
                                if isa(stmt, Expr)
                                    count = 0
                                    for arg in stmt.args
                                        if Base.Meta.isexpr(arg, :replaceglobalref)
                                            count += 1
                                        end
                                    end
                                    count > 0 && return count + 1
                                end
                                return nothing
                            end,
                            (x, i) -> begin
                                items = Any[]
                                stmt = Base.Meta.isexpr(x, :(=)) ? x.args[2] : x
                                if Base.Meta.isexpr(stmt, :replaceglobalref)
                                    tagmodssa = SSAValue(stmt.args[1])
                                    globalref = stmt.args[2]
                                    name = QuoteNode(globalref.name)
                                    result = Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, :tagged_globalref)), overdub_ctx_slot, tagmodssa, name, globalref)
                                elseif isa(stmt, Expr)
                                    result = Expr(stmt.head)
                                    for arg in stmt.args
                                        if Base.Meta.isexpr(arg, :replaceglobalref)
                                            tagmodssa = SSAValue(arg.args[1])
                                            globalref = arg.args[2]
                                            name = QuoteNode(globalref.name)
                                            push!(result.args, SSAValue(i + length(items)))
                                            push!(items, Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, :tagged_globalref)), overdub_ctx_slot, tagmodssa, name, globalref))
                                        else
                                            push!(result.args, arg)
                                        end
                                    end
                                end
                                if Base.Meta.isexpr(x, :(=))
                                    result = Expr(:(=), x.args[1], result)
                                end
                                push!(items, result)
                                return items
                            end)
    end

    #=== untag all `foreigncall` SSAValue/SlotNumber arguments if tagging is enabled ===#

    if istaggingenabled && !iskwfunc
        insert_statements!(overdubbed_code, overdubbed_codelocs,
                            (x, i) -> begin
                                stmt = Base.Meta.isexpr(x, :(=)) ? x.args[2] : x
                                if Base.Meta.isexpr(stmt, :foreigncall)
                                    count = 0
                                    for arg in stmt.args
                                        if isa(arg, SSAValue) || isa(arg, SlotNumber)
                                            count += 1
                                        end
                                    end
                                    count > 0 && return count + 1
                                end
                                return nothing
                            end,
                            (x, i) -> begin
                                items = Any[]
                                stmt = Base.Meta.isexpr(x, :(=)) ? x.args[2] : x
                                result = Expr(:foreigncall)
                                for arg in stmt.args
                                    if isa(arg, SSAValue) || isa(arg, SlotNumber)
                                        push!(result.args, SSAValue(i + length(items)))
                                        push!(items, Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, :untag)), arg, overdub_ctx_slot))
                                    else
                                        push!(result.args, arg)
                                    end
                                end
                                if Base.Meta.isexpr(x, :(=))
                                    result = Expr(:(=), x.args[1], result)
                                end
                                push!(items, result)
                                return items
                            end)
    end

    #=== untag `gotoifnot` conditionals if tagging is enabled ===#

    if istaggingenabled && !iskwfunc
        insert_statements!(overdubbed_code, overdubbed_codelocs,
                            (x, i) -> Base.Meta.isexpr(x, :gotoifnot) ? 2 : nothing,
                            (x, i) -> [
                                Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, :untag)), x.args[1], overdub_ctx_slot),
                                Expr(:gotoifnot, SSAValue(i), x.args[2])
                            ])
    end

    #=== replace `Expr(:new, ...)`/`Expr(:splatnew, ...)` with                            ===#
    #=== `Expr(:call, :tagged_new)`/`Expr(:call, :tagged_splatnew)` if tagging is enabled ===#

    if istaggingenabled && !iskwfunc
        replace_match!(x -> Base.Meta.isexpr(x, :new) || Base.Meta.isexpr(x, :splatnew), overdubbed_code) do x
            tagged_version = x.head == :new ? :tagged_new : :tagged_splatnew
            return Expr(:call, Expr(:nooverdub, GlobalRef(Cassette, tagged_version)), overdub_ctx_slot, x.args...)
        end
    end

    #=== replace `Expr(:call, ...)` with `Expr(:call, :overdub, ...)` calls ===#

    if iskwfunc
        # Another assumption of this `iskwfunc` hack is that the second to last statement in
        # the lowered IR for `Core.kwfunc(f)` is the call to the "underlying" non-kwargs form
        # of `f`. Thus, we `recurse` that call instead of replacing it with `call`.
        for i in 1:length(overdubbed_code)
            stmt = overdubbed_code[i]
            replacewith = i === (length(overdubbed_code) - 1) ? :recurse : :call
            if Base.Meta.isexpr(stmt, :(=))
                replacein = stmt.args
                replaceat = 2
            else
                replacein = overdubbed_code
                replaceat = i
            end
            stmt = replacein[replaceat]
            if Base.Meta.isexpr(stmt, :call)
                replacein[replaceat] = Expr(:call, GlobalRef(Cassette, replacewith), overdub_ctx_slot, stmt.args...)
            end
        end
    else
        arehooksenabled = hashooks(context_type)
        stmtcount = (x, i) -> begin
            i >= original_code_start_index || return nothing
            isassign = Base.Meta.isexpr(x, :(=))
            stmt = isassign ? x.args[2] : x
            if Base.Meta.isexpr(stmt, :call) && !(Base.Meta.isexpr(stmt.args[1], :nooverdub))
                isapplycall = is_ir_element(stmt.args[1], GlobalRef(Core, :_apply), overdubbed_code)
                if isapplycall && arehooksenabled
                    return 7
                elseif isapplycall
                    return 2 + isassign
                elseif arehooksenabled
                    return 4
                else
                    return 1 + isassign
                end
            end
            return nothing
        end
        newstmts = (x, i) -> begin
            callstmt = Base.Meta.isexpr(x, :(=)) ? x.args[2] : x
            isapplycall = is_ir_element(callstmt.args[1], GlobalRef(Core, :_apply), overdubbed_code)
            if isapplycall && arehooksenabled
                callf = callstmt.args[2]
                callargs = callstmt.args[3:end]
                stmts = Any[
                    Expr(:call, GlobalRef(Core, :tuple), overdub_ctx_slot),
                    Expr(:call, GlobalRef(Core, :tuple), callf),
                    Expr(:call, GlobalRef(Core, :_apply), GlobalRef(Cassette, :prehook), SSAValue(i), SSAValue(i + 1), callargs...),
                    Expr(:call, GlobalRef(Core, :_apply), GlobalRef(Cassette, :overdub), SSAValue(i), SSAValue(i + 1), callargs...),
                    Expr(:call, GlobalRef(Core, :tuple), SSAValue(i + 3)),
                    Expr(:call, GlobalRef(Core, :_apply), GlobalRef(Cassette, :posthook), SSAValue(i), SSAValue(i + 4), SSAValue(i + 1), callargs...),
                    Base.Meta.isexpr(x, :(=)) ? Expr(:(=), x.args[1], SSAValue(i + 3)) : SSAValue(i + 3)
                ]
            elseif isapplycall
                callf = callstmt.args[2]
                callargs = callstmt.args[3:end]
                stmts = Any[
                    Expr(:call, GlobalRef(Core, :tuple), overdub_ctx_slot, callf),
                    Expr(:call, GlobalRef(Core, :_apply), GlobalRef(Cassette, :overdub), SSAValue(i), callargs...),
                ]
                Base.Meta.isexpr(x, :(=)) && push!(stmts, Expr(:(=), x.args[1], SSAValue(i + 1)))
            elseif arehooksenabled
                stmts = Any[
                    Expr(:call, GlobalRef(Cassette, :prehook), overdub_ctx_slot, callstmt.args...),
                    Expr(:call, GlobalRef(Cassette, :overdub), overdub_ctx_slot, callstmt.args...),
                    Expr(:call, GlobalRef(Cassette, :posthook), overdub_ctx_slot, SSAValue(i + 1), callstmt.args...),
                    Base.Meta.isexpr(x, :(=)) ? Expr(:(=), x.args[1], SSAValue(i + 1)) : SSAValue(i + 1)
                ]
            else
                stmts = Any[
                    Expr(:call, GlobalRef(Cassette, :overdub), overdub_ctx_slot, callstmt.args...),
                ]
                Base.Meta.isexpr(x, :(=)) && push!(stmts, Expr(:(=), x.args[1], SSAValue(i)))
            end
            return stmts
        end
        insert_statements!(overdubbed_code, overdubbed_codelocs, stmtcount, newstmts)
    end

    #=== unwrap all `Expr(:nooverdub)`s ===#

    replace_match!(x -> x.args[1], x -> Base.Meta.isexpr(x, :nooverdub), overdubbed_code)

    #=== replace all `Expr(:contextslot)`s ===#

    replace_match!(x -> overdub_ctx_slot, x -> Base.Meta.isexpr(x, :contextslot), overdubbed_code)

    #=== set `code_info`/`reflection` fields accordingly ===#

    if code_info.method_for_inference_limit_heuristics === nothing
        code_info.method_for_inference_limit_heuristics = method
    end

    code_info.code = overdubbed_code
    code_info.codelocs = overdubbed_codelocs
    code_info.ssavaluetypes = length(overdubbed_code)
    reflection.code_info = code_info

    return reflection
end

@eval _overdub_fallback($OVERDUB_CONTEXT_NAME, $OVERDUB_ARGUMENTS_NAME...) = fallback($OVERDUB_CONTEXT_NAME, $OVERDUB_ARGUMENTS_NAME...)

const OVERDUB_FALLBACK = begin
    code_info = reflect((typeof(_overdub_fallback), Any, Vararg{Any})).code_info
    code_info.inlineable = true
    code_info
end

# `args` is `(typeof(original_function), map(typeof, original_args_tuple)...)`
function __overdub_generator__(self, context_type, args::Tuple)
    if nfields(args) > 0
        is_builtin = args[1] <: Core.Builtin
        is_invoke = args[1] === typeof(Core.invoke)
        if !is_builtin || is_invoke
            try
                untagged_args = ((untagtype(args[i], context_type) for i in 1:nfields(args))...,)
                reflection = reflect(untagged_args)
                if isa(reflection, Reflection)
                    result = overdub_pass!(reflection, context_type, is_invoke)
                    isa(result, Expr) && return result
                    return reflection.code_info
                end
            catch err
                errmsg = "ERROR COMPILING $args IN CONTEXT $(context_type): \n" * sprint(showerror, err)
                errmsg *= "\n" .* repr("text/plain", stacktrace(catch_backtrace()))
                return quote
                    error($errmsg)
                end
            end
        end
    end
    return copy_code_info(OVERDUB_FALLBACK)
end

function overdub end

function recurse end

recurse(ctx::Context, ::typeof(Core._apply), f, args...) = Core._apply(recurse, (ctx, f), args...)

function delete_overdub_or_recurse_method(overdub_or_recurse)
    ms = Base.methods(overdub_or_recurse, (Context, Any)).ms
    for m in ms
        if m.sig == Tuple{typeof(overdub_or_recurse),Context,Vararg{Any}}
            Base.delete_method(m)
            return nothing
        end
    end
    return nothing
end

function overdub_definition(line, file)
    return quote
        $Cassette.delete_overdub_or_recurse_method($Cassette.overdub)
        $Cassette.delete_overdub_or_recurse_method($Cassette.recurse)
        function $Cassette.overdub($OVERDUB_CONTEXT_NAME::$Cassette.Context, $OVERDUB_ARGUMENTS_NAME...)
            $(Expr(:meta, :generated_only))
            $(Expr(:meta,
                   :generated,
                   Expr(:new,
                        Core.GeneratedFunctionStub,
                        :__overdub_generator__,
                        Any[:overdub, OVERDUB_CONTEXT_NAME, OVERDUB_ARGUMENTS_NAME],
                        Any[],
                        line,
                        QuoteNode(Symbol(file)),
                        true)))
        end
        function $Cassette.recurse($OVERDUB_CONTEXT_NAME::$Cassette.Context, $OVERDUB_ARGUMENTS_NAME...)
            $(Expr(:meta, :generated_only))
            $(Expr(:meta,
                   :generated,
                   Expr(:new,
                        Core.GeneratedFunctionStub,
                        :__overdub_generator__,
                        Any[:recurse, OVERDUB_CONTEXT_NAME, OVERDUB_ARGUMENTS_NAME],
                        Any[],
                        line,
                        QuoteNode(Symbol(file)),
                        true)))
        end
    end
end

@doc(
"""
```
overdub(context::Context, f, args...)
```

Execute `f(args...)` overdubbed with respect to `context`.

More specifically, execute `f(args...)`, but with every internal method invocation `g(x...)`
replaced by statements similar to the following:

```
begin
    prehook(context, g, x...)
    overdub(context, g, x...) # %n
    posthook(context, %n, g, x...)
    %n
end
```

Otherwise, if Cassette cannot retrieve lowered IR for the method body of `f(args...)`,
then `fallback(context, f, args...)` will be called instead. Cassette's [`canrecurse`](@ref)
function is a useful utility for checking if this will occur.

If the injected `prehook`/`posthook` statements are not needed for your use
case, you can disable their injection via the [`disablehooks`](@ref) function.

Additionally, for every method body encountered in the execution trace, apply the compiler pass
associated with `context` if one exists. Note that this user-provided pass is performed on
the method IR before method invocations are transformed into the form specified above. See
the [`@pass`](@ref) macro for further details.

If `Cassette.hastagging(typeof(context))`, then a number of additional passes are run in
order to accomodate tagged value propagation:

- `Expr(:new)` is replaced with a call to `Cassette.tagged_new`
- `Expr(:splatnew)` is replaced with a call to `Cassette.tagged_splatnew`
- conditional values passed to `Expr(:gotoifnot)` are untagged
- arguments to `Expr(:foreigncall)` are untagged
- load/stores to external module bindings are intercepted by the tagging system

The default definition of `overdub` is to recursively enter the given function
and continue overdubbing, but one can interrupt/redirect this recursion by
overloading `overdub` w.r.t. a given context and/or method signature to define
new contextual execution primitives. For example:

```
julia> using Cassette

julia> Cassette.@context Ctx;

julia> Cassette.overdub(::Ctx, ::typeof(sin), x) = cos(x)

julia> Cassette.overdub(Ctx(), x -> sin(x) + cos(x), 1) == 2 * cos(1)
true
```

See also: [`recurse`](@ref), [`prehook`](@ref), [`posthook`](@ref)
""",
overdub)

@doc(
"""
```
recurse(context::Context, f, args...)
```

Execute `f(args...)` overdubbed with respect to `context`.

This method performs exactly the same transformation as the default
[`overdub`](@ref) transformation, but is not meant to be overloaded. Thus, one
can call `recurse` to "continue" recursively overdubbing a function when calling
`overdub` directly on that function might've dispatched to a contextual
primitive.

To illustrate why `recurse` might be useful, consider the following example
which utilizes `recurse` as part of a Cassette-based memoization implementation
for the classic Fibonacci function:

```
using Cassette: Cassette, @context, overdub, recurse

fib(x) = x < 3 ? 1 : fib(x - 2) + fib(x - 1)
fibtest(n) = fib(2 * n) + n

@context MemoizeCtx

function Cassette.overdub(ctx::MemoizeCtx, ::typeof(fib), x)
    result = get(ctx.metadata, x, 0)
    if result === 0
        result = recurse(ctx, fib, x)
        ctx.metadata[x] = result
    end
    return result
end
```

See Cassette's Contextual Dispatch documentation for more details and examples.
""",
recurse)

"""
```
Cassette.@overdub(ctx, expression)
```

A convenience macro for executing `expression` within the context `ctx`. This
macro roughly expands to `Cassette.recurse(ctx, () -> expression)`.

See also: [`overdub`](@ref), [`recurse`](@ref)
"""
macro overdub(ctx, expr)
    return :(recurse($(esc(ctx)), () -> $(esc(expr))))
end


# This file is a part of Julia. License is MIT: https://julialang.org/license

# Operations with the file system (paths) ##

export
    cd,
    chmod,
    chown,
    cp,
    cptree,
    mkdir,
    mkpath,
    mktemp,
    mktempdir,
    mv,
    pwd,
    rename,
    readlink,
    readdir,
    rm,
    samefile,
    sendfile,
    symlink,
    tempdir,
    tempname,
    touch,
    unlink,
    walkdir

# get and set current directory

"""
    pwd() -> AbstractString

Get the current working directory.

# Examples
```julia-repl
julia> pwd()
"/home/JuliaUser"

julia> cd("/home/JuliaUser/Projects/julia")

julia> pwd()
"/home/JuliaUser/Projects/julia"
```
"""
function pwd()
    buf = Base.StringVector(AVG_PATH - 1) # space for null-terminator implied by StringVector
    sz = RefValue{Csize_t}(length(buf) + 1) # total buffer size including null
    while true
        rc = ccall(:uv_cwd, Cint, (Ptr{UInt8}, Ptr{Csize_t}), buf, sz)
        if rc == 0
            resize!(buf, sz[])
            return String(buf)
        elseif rc == Base.UV_ENOBUFS
            resize!(buf, sz[] - 1) # space for null-terminator implied by StringVector
        else
            uv_error(:cwd, rc)
        end
    end
end


"""
    cd(dir::AbstractString=homedir())

Set the current working directory.

# Examples
```julia-repl
julia> cd("/home/JuliaUser/Projects/julia")

julia> pwd()
"/home/JuliaUser/Projects/julia"

julia> cd()

julia> pwd()
"/home/JuliaUser"
```
"""
function cd(dir::AbstractString)
    uv_error("chdir $dir", ccall(:uv_chdir, Cint, (Cstring,), dir))
end
cd() = cd(homedir())

if Sys.iswindows()
    function cd(f::Function, dir::AbstractString)
        old = pwd()
        try
            cd(dir)
            f()
       finally
            cd(old)
        end
    end
else
    function cd(f::Function, dir::AbstractString)
        fd = ccall(:open, Int32, (Cstring, Int32), :., 0)
        systemerror(:open, fd == -1)
        try
            cd(dir)
            f()
        finally
            systemerror(:fchdir, ccall(:fchdir, Int32, (Int32,), fd) != 0)
            systemerror(:close, ccall(:close, Int32, (Int32,), fd) != 0)
        end
    end
end
"""
    cd(f::Function, dir::AbstractString=homedir())

Temporarily change the current working directory to `dir`, apply function `f` and
finally return to the original directory.

# Examples
```julia-repl
julia> pwd()
"/home/JuliaUser"

julia> cd(readdir, "/home/JuliaUser/Projects/julia")
34-element Array{String,1}:
 ".circleci"
 ".freebsdci.sh"
 ".git"
 ".gitattributes"
 ".github"
 ⋮
 "test"
 "ui"
 "usr"
 "usr-staging"

julia> pwd()
"/home/JuliaUser"
```
"""
cd(f::Function) = cd(f, homedir())

function checkmode(mode::Integer)
    if !(0 <= mode <= 511)
        throw(ArgumentError("Mode must be between 0 and 511 = 0o777"))
    end
    mode
end

"""
    mkdir(path::AbstractString; mode::Unsigned = 0o777)

Make a new directory with name `path` and permissions `mode`. `mode` defaults to `0o777`,
modified by the current file creation mask. This function never creates more than one
directory. If the directory already exists, or some intermediate directories do not exist,
this function throws an error. See [`mkpath`](@ref) for a function which creates all
required intermediate directories.
Return `path`.

# Examples
```julia-repl
julia> mkdir("testingdir")
"testingdir"

julia> cd("testingdir")

julia> pwd()
"/home/JuliaUser/testingdir"
```
"""
function mkdir(path::AbstractString; mode::Integer = 0o777)
    req = Libc.malloc(_sizeof_uv_fs)
    try
        ret = ccall(:uv_fs_mkdir, Cint,
                    (Ptr{Cvoid}, Ptr{Cvoid}, Cstring, Cint, Ptr{Cvoid}),
                    C_NULL, req, path, checkmode(mode), C_NULL)
        if ret < 0
            ccall(:uv_fs_req_cleanup, Cvoid, (Ptr{Cvoid},), req)
            uv_error("mkdir", ret)
        end
        ccall(:uv_fs_req_cleanup, Cvoid, (Ptr{Cvoid},), req)
        return path
    finally
        Libc.free(req)
    end
end

"""
    mkpath(path::AbstractString; mode::Unsigned = 0o777)

Create all directories in the given `path`, with permissions `mode`. `mode` defaults to
`0o777`, modified by the current file creation mask.
Return `path`.

# Examples
```julia-repl
julia> mkdir("testingdir")
"testingdir"

julia> cd("testingdir")

julia> pwd()
"/home/JuliaUser/testingdir"

julia> mkpath("my/test/dir")
"my/test/dir"

julia> readdir()
1-element Array{String,1}:
 "my"

julia> cd("my")

julia> readdir()
1-element Array{String,1}:
 "test"

julia> readdir("test")
1-element Array{String,1}:
 "dir"
```
"""
function mkpath(path::AbstractString; mode::Integer = 0o777)
    isdirpath(path) && (path = dirname(path))
    dir = dirname(path)
    (path == dir || isdir(path)) && return path
    mkpath(dir, mode = checkmode(mode))
    try
        mkdir(path, mode = mode)
    catch err
        # If there is a problem with making the directory, but the directory
        # does in fact exist, then ignore the error. Else re-throw it.
        if !isa(err, SystemError) || !isdir(path)
            rethrow()
        end
    end
    path
end

"""
    rm(path::AbstractString; force::Bool=false, recursive::Bool=false)

Delete the file, link, or empty directory at the given path. If `force=true` is passed, a
non-existing path is not treated as error. If `recursive=true` is passed and the path is a
directory, then all contents are removed recursively.

# Examples
```jldoctest
julia> mkpath("my/test/dir");

julia> rm("my", recursive=true)

julia> rm("this_file_does_not_exist", force=true)

julia> rm("this_file_does_not_exist")
ERROR: IOError: unlink: no such file or directory (ENOENT)
Stacktrace:
[...]
```
"""
function rm(path::AbstractString; force::Bool=false, recursive::Bool=false)
    if islink(path) || !isdir(path)
        try
            @static if Sys.iswindows()
                # is writable on windows actually means "is deletable"
                if (filemode(lstat(path)) & 0o222) == 0
                    chmod(path, 0o777)
                end
            end
            unlink(path)
        catch err
            if force && isa(err, IOError) && err.code==Base.UV_ENOENT
                return
            end
            rethrow()
        end
    else
        if recursive
            for p in readdir(path)
                rm(joinpath(path, p), force=force, recursive=true)
            end
        end
        @static if Sys.iswindows()
            ret = ccall(:_wrmdir, Int32, (Cwstring,), path)
        else
            ret = ccall(:rmdir, Int32, (Cstring,), path)
        end
        systemerror(:rmdir, ret != 0, extrainfo=path)
    end
end


# The following use Unix command line facilities
function checkfor_mv_cp_cptree(src::AbstractString, dst::AbstractString, txt::AbstractString;
                                                          force::Bool=false)
    if ispath(dst)
        if force
            # Check for issue when: (src == dst) or when one is a link to the other
            # https://github.com/JuliaLang/julia/pull/11172#issuecomment-100391076
            if Base.samefile(src, dst)
                abs_src = islink(src) ? abspath(readlink(src)) : abspath(src)
                abs_dst = islink(dst) ? abspath(readlink(dst)) : abspath(dst)
                throw(ArgumentError(string("'src' and 'dst' refer to the same file/dir.",
                                           "This is not supported.\n  ",
                                           "`src` refers to: $(abs_src)\n  ",
                                           "`dst` refers to: $(abs_dst)\n")))
            end
            rm(dst; recursive=true)
        else
            throw(ArgumentError(string("'$dst' exists. `force=true` ",
                                       "is required to remove '$dst' before $(txt).")))
        end
    end
end

function cptree(src::AbstractString, dst::AbstractString; force::Bool=false,
                                                          follow_symlinks::Bool=false)
    isdir(src) || throw(ArgumentError("'$src' is not a directory. Use `cp(src, dst)`"))
    checkfor_mv_cp_cptree(src, dst, "copying"; force=force)
    mkdir(dst)
    for name in readdir(src)
        srcname = joinpath(src, name)
        if !follow_symlinks && islink(srcname)
            symlink(readlink(srcname), joinpath(dst, name))
        elseif isdir(srcname)
            cptree(srcname, joinpath(dst, name); force=force,
                                                 follow_symlinks=follow_symlinks)
        else
            sendfile(srcname, joinpath(dst, name))
        end
    end
end

"""
    cp(src::AbstractString, dst::AbstractString; force::Bool=false, follow_symlinks::Bool=false)

Copy the file, link, or directory from `src` to `dst`.
`force=true` will first remove an existing `dst`.

If `follow_symlinks=false`, and `src` is a symbolic link, `dst` will be created as a
symbolic link. If `follow_symlinks=true` and `src` is a symbolic link, `dst` will be a copy
of the file or directory `src` refers to.
Return `dst`.
"""
function cp(src::AbstractString, dst::AbstractString; force::Bool=false,
                                                      follow_symlinks::Bool=false)
    checkfor_mv_cp_cptree(src, dst, "copying"; force=force)
    if !follow_symlinks && islink(src)
        symlink(readlink(src), dst)
    elseif isdir(src)
        cptree(src, dst; force=force, follow_symlinks=follow_symlinks)
    else
        sendfile(src, dst)
    end
    dst
end

"""
    mv(src::AbstractString, dst::AbstractString; force::Bool=false)

Move the file, link, or directory from `src` to `dst`.
`force=true` will first remove an existing `dst`.
Return `dst`.

# Examples
```jldoctest; filter = r"Stacktrace:(\\n \\[[0-9]+\\].*)*"
julia> write("hello.txt", "world");

julia> mv("hello.txt", "goodbye.txt")
"goodbye.txt"

julia> "hello.txt" in readdir()
false

julia> readline("goodbye.txt")
"world"

julia> write("hello.txt", "world2");

julia> mv("hello.txt", "goodbye.txt")
ERROR: ArgumentError: 'goodbye.txt' exists. `force=true` is required to remove 'goodbye.txt' before moving.
Stacktrace:
 [1] #checkfor_mv_cp_cptree#10(::Bool, ::Function, ::String, ::String, ::String) at ./file.jl:293
[...]

julia> mv("hello.txt", "goodbye.txt", force=true)
"goodbye.txt"

julia> rm("goodbye.txt");

```
"""
function mv(src::AbstractString, dst::AbstractString; force::Bool=false)
    checkfor_mv_cp_cptree(src, dst, "moving"; force=force)
    rename(src, dst)
    dst
end

"""
    touch(path::AbstractString)

Update the last-modified timestamp on a file to the current time.

If the file does not exist a new file is created.

Return `path`.

# Examples
```julia-repl
julia> write("my_little_file", 2);

julia> mtime("my_little_file")
1.5273815391135583e9

julia> touch("my_little_file");

julia> mtime("my_little_file")
1.527381559163435e9
```

We can see the [`mtime`](@ref) has been modified by `touch`.
"""
function touch(path::AbstractString)
    f = open(path, JL_O_WRONLY | JL_O_CREAT, 0o0666)
    try
        if Sys.isunix()
            ret = ccall(:futimes, Cint, (Cint, Ptr{Cvoid}), fd(f), C_NULL)
            systemerror(:futimes, ret != 0, extrainfo=path)
        else
            t = time()
            futime(f,t,t)
        end
    finally
        close(f)
    end
    path
end

"""
    tempdir()

Gets the path of the temporary directory. On Windows, `tempdir()` uses the first environment
variable found in the ordered list `TMP`, `TEMP`, `USERPROFILE`. On all other operating
systems, `tempdir()` uses the first environment variable found in the ordered list `TMPDIR`,
`TMP`, `TEMP`, and `TEMPDIR`. If none of these are found, the path `"/tmp"` is used.
"""
function tempdir()
    buf = Base.StringVector(AVG_PATH - 1) # space for null-terminator implied by StringVector
    sz = RefValue{Csize_t}(length(buf) + 1) # total buffer size including null
    while true
        rc = ccall(:uv_os_tmpdir, Cint, (Ptr{UInt8}, Ptr{Csize_t}), buf, sz)
        if rc == 0
            resize!(buf, sz[])
            return String(buf)
        elseif rc == Base.UV_ENOBUFS
            resize!(buf, sz[] - 1)  # space for null-terminator implied by StringVector
        else
            uv_error(:tmpdir, rc)
        end
    end
end

const TEMP_CLEANUP_MIN = Ref(1024)
const TEMP_CLEANUP_MAX = Ref(1024)
const TEMP_CLEANUP = Dict{String,Bool}()
const TEMP_CLEANUP_LOCK = ReentrantLock()

function temp_cleanup_later(path::AbstractString; asap::Bool=false)
    lock(TEMP_CLEANUP_LOCK)
    # each path should only be inserted here once, but if there
    # is a collision, let !asap win over asap: if any user might
    # still be using the path, don't delete it until process exit
    TEMP_CLEANUP[path] = get(TEMP_CLEANUP, path, true) & asap
    if length(TEMP_CLEANUP) > TEMP_CLEANUP_MAX[]
        temp_cleanup_purge(false)
        TEMP_CLEANUP_MAX[] = max(TEMP_CLEANUP_MIN[], 2*length(TEMP_CLEANUP))
    end
    unlock(TEMP_CLEANUP_LOCK)
    return nothing
end

function temp_cleanup_purge(all::Bool=true)
    need_gc = Sys.iswindows()
    for (path, asap) in TEMP_CLEANUP
        try
            if (all || asap) && ispath(path)
                need_gc && GC.gc(true)
                need_gc = false
                rm(path, recursive=true, force=true)
            end
            !ispath(path) && delete!(TEMP_CLEANUP, path)
        catch ex
            @warn "temp cleanup" _group=:file exception=(ex, catch_backtrace())
        end
    end
end

const temp_prefix = "jl_"

if Sys.iswindows()

function _win_tempname(temppath::AbstractString, uunique::UInt32)
    tempp = cwstring(temppath)
    temppfx = cwstring(temp_prefix)
    tname = Vector{UInt16}(undef, 32767)
    uunique = ccall(:GetTempFileNameW, stdcall, UInt32,
                    (Ptr{UInt16}, Ptr{UInt16}, UInt32, Ptr{UInt16}),
                    tempp, temppfx, uunique, tname)
    windowserror("GetTempFileName", uunique == 0)
    lentname = something(findfirst(iszero, tname))
    @assert lentname > 0
    resize!(tname, lentname - 1)
    return transcode(String, tname)
end

function mktemp(parent::AbstractString=tempdir(); cleanup::Bool=true)
    filename = _win_tempname(parent, UInt32(0))
    cleanup && temp_cleanup_later(filename)
    return (filename, Base.open(filename, "r+"))
end

# generate a random string from random bytes
function _rand_string()
    nchars = 10
    A = Vector{UInt8}(undef, nchars)
    windowserror("SystemFunction036 (RtlGenRandom)", 0 == ccall(
        (:SystemFunction036, :Advapi32), stdcall, UInt8, (Ptr{Cvoid}, UInt32),
            A, sizeof(A)))

    slug = Base.StringVector(10)
    chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i = 1:nchars
        slug[i] = chars[(A[i] % length(chars)) + 1]
    end
    return name = String(slug)
end

function tempname(parent::AbstractString=tempdir(); cleanup::Bool=true)
    isdir(parent) || throw(ArgumentError("$(repr(parent)) is not a directory"))
    name = _rand_string()
    filename = joinpath(parent, temp_prefix * name)
    @assert !ispath(filename)
    cleanup && temp_cleanup_later(filename)
    return filename
end

else # !windows

# Obtain a temporary filename.
function tempname(parent::AbstractString=tempdir(); cleanup::Bool=true)
    isdir(parent) || throw(ArgumentError("$(repr(parent)) is not a directory"))
    p = ccall(:tempnam, Cstring, (Cstring, Cstring), parent, temp_prefix)
    systemerror(:tempnam, p == C_NULL)
    s = unsafe_string(p)
    Libc.free(p)
    cleanup && temp_cleanup_later(s)
    return s
end

# Create and return the name of a temporary file along with an IOStream
function mktemp(parent::AbstractString=tempdir(); cleanup::Bool=true)
    b = joinpath(parent, temp_prefix * "XXXXXX")
    p = ccall(:mkstemp, Int32, (Cstring,), b) # modifies b
    systemerror(:mktemp, p == -1)
    cleanup && temp_cleanup_later(b)
    return (b, fdio(p, true))
end


end # os-test


"""
    tempname(parent=tempdir(); cleanup=true) -> String

Generate a temporary file path. This function only returns a path; no file is
created. The path is likely to be unique, but this cannot be guaranteed due to
the very remote posibility of two simultaneous calls to `tempname` generating
the same file name. The name is guaranteed to differ from all files already
existing at the time of the call to `tempname`.

When called with no arguments, the temporary name will be an absolute path to a
temporary name in the system temporary directory as given by `tempdir()`. If a
`parent` directory argument is given, the temporary path will be in that
directory instead.

The `cleanup` option controls whether the process attempts to delete the
returned path automatically when the process exits. Note that the `tempname`
function does not create any file or directory at the returned location, so
there is nothing to cleanup unless you create a file or directory there. If
you do and `clean` is `true` it will be deleted upon process termination.

!!! compat "Julia 1.4"
    The `parent` and `cleanup` arguments were added in 1.4. Prior to Julia 1.4
    the path `tempname` would never be cleaned up at process termination.

!!! warning

    This can lead to security holes if another process obtains the same
    file name and creates the file before you are able to. Open the file with
    `JL_O_EXCL` if this is a concern. Using [`mktemp()`](@ref) is also
    recommended instead.
"""
tempname()

"""
    mktemp(parent=tempdir(); cleanup=true) -> (path, io)

Return `(path, io)`, where `path` is the path of a new temporary file in `parent`
and `io` is an open file object for this path. The `cleanup` option controls whether
the temporary file is automatically deleted when the process exits.
"""
mktemp(parent)

"""
    mktempdir(parent=tempdir(); prefix=$(repr(temp_prefix)), cleanup=true) -> path

Create a temporary directory in the `parent` directory with a name
constructed from the given prefix and a random suffix, and return its path.
Additionally, any trailing `X` characters may be replaced with random characters.
If `parent` does not exist, throw an error. The `cleanup` option controls whether
the temporary directory is automatically deleted when the process exits.
"""
function mktempdir(parent::AbstractString=tempdir();
    prefix::AbstractString=temp_prefix, cleanup::Bool=true)
    if isempty(parent) || occursin(path_separator_re, parent[end:end])
        # append a path_separator only if parent didn't already have one
        tpath = "$(parent)$(prefix)XXXXXX"
    else
        tpath = "$(parent)$(path_separator)$(prefix)XXXXXX"
    end

    req = Libc.malloc(_sizeof_uv_fs)
    try
        ret = ccall(:uv_fs_mkdtemp, Cint,
                    (Ptr{Cvoid}, Ptr{Cvoid}, Cstring, Ptr{Cvoid}),
                    C_NULL, req, tpath, C_NULL)
        if ret < 0
            ccall(:uv_fs_req_cleanup, Cvoid, (Ptr{Cvoid},), req)
            uv_error("mktempdir", ret)
        end
        path = unsafe_string(ccall(:jl_uv_fs_t_path, Cstring, (Ptr{Cvoid},), req))
        ccall(:uv_fs_req_cleanup, Cvoid, (Ptr{Cvoid},), req)
        cleanup && temp_cleanup_later(path)
        return path
    finally
        Libc.free(req)
    end
end


"""
    mktemp(f::Function, parent=tempdir())

Apply the function `f` to the result of [`mktemp(parent)`](@ref) and remove the
temporary file upon completion.
"""
function mktemp(fn::Function, parent::AbstractString=tempdir())
    (tmp_path, tmp_io) = mktemp(parent, cleanup=false)
    try
        fn(tmp_path, tmp_io)
    finally
        try
            close(tmp_io)
            ispath(tmp_path) && rm(tmp_path)
        catch ex
            @error "mktemp cleanup" _group=:file exception=(ex, catch_backtrace())
            # might be possible to remove later
            temp_cleanup_later(tmp_path, asap=true)
        end
    end
end

"""
    mktempdir(f::Function, parent=tempdir(); prefix=$(repr(temp_prefix)))

Apply the function `f` to the result of [`mktempdir(parent; prefix)`](@ref) and remove the
temporary directory all of its contents upon completion.
"""
function mktempdir(fn::Function, parent::AbstractString=tempdir();
    prefix::AbstractString=temp_prefix)
    tmpdir = mktempdir(parent; prefix=prefix, cleanup=false)
    try
        fn(tmpdir)
    finally
        try
            ispath(tmpdir) && rm(tmpdir, recursive=true)
        catch ex
            @error "mktempdir cleanup" _group=:file exception=(ex, catch_backtrace())
            # might be possible to remove later
            temp_cleanup_later(tmpdir, asap=true)
        end
    end
end

struct uv_dirent_t
    name::Ptr{UInt8}
    typ::Cint
end

"""
    readdir(dir::AbstractString=pwd();
        join::Bool = false,
        sort::Bool = true,
    ) -> Vector{String}

Return the names in the directory `dir` or the current working directory if not
given. When `join` is false, `readdir` returns just the names in the directory
as is; when `join` is true, it returns `joinpath(dir, name)` for each `name` so
that the returned strings are full paths. If you want to get absolute paths
back, call `readdir` with an absolute directory path and `join` set to true.

By default, `readdir` sorts the list of names it returns. If you want to skip
sorting the names and get them in the order that the file system lists them,
you can use `readir(dir, sort=false)` to opt out of sorting.

!!! compat "Julia 1.4"
    The `join` and `sort` keyword arguments require at least Julia 1.4.

# Examples
```julia-repl
julia> cd("/home/JuliaUser/dev/julia")

julia> readdir()
30-element Array{String,1}:
 ".appveyor.yml"
 ".git"
 ".gitattributes"
 ⋮
 "ui"
 "usr"
 "usr-staging"

julia> readdir(join=true)
30-element Array{String,1}:
 "/home/JuliaUser/dev/julia/.appveyor.yml"
 "/home/JuliaUser/dev/julia/.git"
 "/home/JuliaUser/dev/julia/.gitattributes"
 ⋮
 "/home/JuliaUser/dev/julia/ui"
 "/home/JuliaUser/dev/julia/usr"
 "/home/JuliaUser/dev/julia/usr-staging"

julia> readdir("base")
145-element Array{String,1}:
 ".gitignore"
 "Base.jl"
 "Enums.jl"
 ⋮
 "version_git.sh"
 "views.jl"
 "weakkeydict.jl"

julia> readdir("base", join=true)
145-element Array{String,1}:
 "base/.gitignore"
 "base/Base.jl"
 "base/Enums.jl"
 ⋮
 "base/version_git.sh"
 "base/views.jl"
 "base/weakkeydict.jl"```

julia> readdir(abspath("base"), join=true)
145-element Array{String,1}:
 "/home/JuliaUser/dev/julia/base/.gitignore"
 "/home/JuliaUser/dev/julia/base/Base.jl"
 "/home/JuliaUser/dev/julia/base/Enums.jl"
 ⋮
 "/home/JuliaUser/dev/julia/base/version_git.sh"
 "/home/JuliaUser/dev/julia/base/views.jl"
 "/home/JuliaUser/dev/julia/base/weakkeydict.jl"
```
"""
function readdir(dir::AbstractString; join::Bool=false, sort::Bool=true)
    # Allocate space for uv_fs_t struct
    uv_readdir_req = zeros(UInt8, ccall(:jl_sizeof_uv_fs_t, Int32, ()))

    # defined in sys.c, to call uv_fs_readdir, which sets errno on error.
    err = ccall(:uv_fs_scandir, Int32, (Ptr{Cvoid}, Ptr{UInt8}, Cstring, Cint, Ptr{Cvoid}),
                C_NULL, uv_readdir_req, dir, 0, C_NULL)
    err < 0 && throw(SystemError("unable to read directory $dir", -err))

    # iterate the listing into entries
    entries = String[]
    ent = Ref{uv_dirent_t}()
    while Base.UV_EOF != ccall(:uv_fs_scandir_next, Cint, (Ptr{Cvoid}, Ptr{uv_dirent_t}), uv_readdir_req, ent)
        name = unsafe_string(ent[].name)
        push!(entries, join ? joinpath(dir, name) : name)
    end

    # Clean up the request string
    ccall(:uv_fs_req_cleanup, Cvoid, (Ptr{UInt8},), uv_readdir_req)

    # sort entries unless opted out
    sort && sort!(entries)

    return entries
end
readdir(; join::Bool=false, sort::Bool=true) =
    readdir(join ? pwd() : ".", join=join, sort=sort)

"""
    walkdir(dir; topdown=true, follow_symlinks=false, onerror=throw)

Return an iterator that walks the directory tree of a directory.
The iterator returns a tuple containing `(rootpath, dirs, files)`.
The directory tree can be traversed top-down or bottom-up.
If `walkdir` encounters a [`SystemError`](@ref)
it will rethrow the error by default.
A custom error handling function can be provided through `onerror` keyword argument.
`onerror` is called with a `SystemError` as argument.

# Examples
```julia
for (root, dirs, files) in walkdir(".")
    println("Directories in \$root")
    for dir in dirs
        println(joinpath(root, dir)) # path to directories
    end
    println("Files in \$root")
    for file in files
        println(joinpath(root, file)) # path to files
    end
end
```

```julia-repl
julia> mkpath("my/test/dir");

julia> itr = walkdir("my");

julia> (root, dirs, files) = first(itr)
("my", ["test"], String[])

julia> (root, dirs, files) = first(itr)
("my/test", ["dir"], String[])

julia> (root, dirs, files) = first(itr)
("my/test/dir", String[], String[])
```
"""
function walkdir(root; topdown=true, follow_symlinks=false, onerror=throw)
    content = nothing
    try
        content = readdir(root)
    catch err
        isa(err, SystemError) || throw(err)
        onerror(err)
        # Need to return an empty closed channel to skip the current root folder
        chnl = Channel(0)
        close(chnl)
        return chnl
    end
    dirs = Vector{eltype(content)}()
    files = Vector{eltype(content)}()
    for name in content
        path = joinpath(root, name)

        # If we're not following symlinks, then treat all symlinks as files
        if (!follow_symlinks && islink(path)) || !isdir(path)
            push!(files, name)
        else
            push!(dirs, name)
        end
    end

    function _it(chnl)
        if topdown
            put!(chnl, (root, dirs, files))
        end
        for dir in dirs
            path = joinpath(root,dir)
            if follow_symlinks || !islink(path)
                for (root_l, dirs_l, files_l) in walkdir(path, topdown=topdown, follow_symlinks=follow_symlinks, onerror=onerror)
                    put!(chnl, (root_l, dirs_l, files_l))
                end
            end
        end
        if !topdown
            put!(chnl, (root, dirs, files))
        end
    end

    return Channel(_it)
end

function unlink(p::AbstractString)
    err = ccall(:jl_fs_unlink, Int32, (Cstring,), p)
    uv_error("unlink", err)
    nothing
end

# For move command
function rename(src::AbstractString, dst::AbstractString)
    err = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), src, dst)
    # on error, default to cp && rm
    if err < 0
        # force: is already done in the mv function
        cp(src, dst; force=false, follow_symlinks=false)
        rm(src; recursive=true)
    end
    nothing
end

function sendfile(src::AbstractString, dst::AbstractString)
    src_open = false
    dst_open = false
    local src_file, dst_file
    try
        src_file = open(src, JL_O_RDONLY)
        src_open = true
        dst_file = open(dst, JL_O_CREAT | JL_O_TRUNC | JL_O_WRONLY, filemode(src_file))
        dst_open = true

        bytes = filesize(stat(src_file))
        sendfile(dst_file, src_file, Int64(0), Int(bytes))
    finally
        if src_open && isopen(src_file)
            close(src_file)
        end
        if dst_open && isopen(dst_file)
            close(dst_file)
        end
    end
end

if Sys.iswindows()
    const UV_FS_SYMLINK_JUNCTION = 0x0002
end

"""
    symlink(target::AbstractString, link::AbstractString)

Creates a symbolic link to `target` with the name `link`.

!!! note
    This function raises an error under operating systems that do not support
    soft symbolic links, such as Windows XP.
"""
function symlink(p::AbstractString, np::AbstractString)
    @static if Sys.iswindows()
        if Sys.windows_version() < Sys.WINDOWS_VISTA_VER
            error("Windows XP does not support soft symlinks")
        end
    end
    flags = 0
    @static if Sys.iswindows()
        if isdir(p)
            flags |= UV_FS_SYMLINK_JUNCTION
            p = abspath(p)
        end
    end
    err = ccall(:jl_fs_symlink, Int32, (Cstring, Cstring, Cint), p, np, flags)
    @static if Sys.iswindows()
        if err < 0 && !isdir(p)
            @warn "On Windows, creating file symlinks requires Administrator privileges" maxlog=1 _group=:file
        end
    end
    uv_error("symlink",err)
end

"""
    readlink(path::AbstractString) -> AbstractString

Return the target location a symbolic link `path` points to.
"""
function readlink(path::AbstractString)
    req = Libc.malloc(_sizeof_uv_fs)
    try
        ret = ccall(:uv_fs_readlink, Int32,
            (Ptr{Cvoid}, Ptr{Cvoid}, Cstring, Ptr{Cvoid}),
            C_NULL, req, path, C_NULL)
        if ret < 0
            ccall(:uv_fs_req_cleanup, Cvoid, (Ptr{Cvoid},), req)
            uv_error("readlink", ret)
            @assert false
        end
        tgt = unsafe_string(ccall(:jl_uv_fs_t_ptr, Cstring, (Ptr{Cvoid},), req))
        ccall(:uv_fs_req_cleanup, Cvoid, (Ptr{Cvoid},), req)
        return tgt
    finally
        Libc.free(req)
    end
end

"""
    chmod(path::AbstractString, mode::Integer; recursive::Bool=false)

Change the permissions mode of `path` to `mode`. Only integer `mode`s (e.g. `0o777`) are
currently supported. If `recursive=true` and the path is a directory all permissions in
that directory will be recursively changed.
Return `path`.
"""
function chmod(path::AbstractString, mode::Integer; recursive::Bool=false)
    err = ccall(:jl_fs_chmod, Int32, (Cstring, Cint), path, mode)
    uv_error("chmod", err)
    if recursive && isdir(path)
        for p in readdir(path)
            if !islink(joinpath(path, p))
                chmod(joinpath(path, p), mode, recursive=true)
            end
        end
    end
    path
end

"""
    chown(path::AbstractString, owner::Integer, group::Integer=-1)

Change the owner and/or group of `path` to `owner` and/or `group`. If the value entered for `owner` or `group`
is `-1` the corresponding ID will not change. Only integer `owner`s and `group`s are currently supported.
Return `path`.
"""
function chown(path::AbstractString, owner::Integer, group::Integer=-1)
    err = ccall(:jl_fs_chown, Int32, (Cstring, Cint, Cint), path, owner, group)
    uv_error("chown",err)
    path
end


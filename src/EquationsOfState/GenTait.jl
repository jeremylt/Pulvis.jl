# ------------------------------------------------------------------------------
# Equations of state - gen tait
# ------------------------------------------------------------------------------

"""
```julia
densitytopressure_GenTait(ρ, ρ0, c, refdensityfactor, exponent)
```

Convert density to pressure

# Arguments:

  - `ρ::Float`:                particle density
  - `ρ0::Float`:               initial density
  - `c::Float`:                speed of sound
  - `refdensityfactor::Float`: reference density factor
  - `exponent::Int`:           exponent

# Returns:

  - particle pressure

# Example:

```jldoctest
# set equation of state
setequationsofstate(EquationOfState.GenTait)

# setup
c = 3.5
refdensityfactor = 0.9
ρ = 0.75
ρ0 = 0.4

# verify
exponent = 7
@assert densitytopressure(ρ, ρ0, c, refdensityfactor, exponent) ≈ 0.5640046918e2
exponent = 1
@assert densitytopressure(ρ, ρ0, c, refdensityfactor, exponent) ≈ 0.47775e1

# output

```
"""
function densitytopressure_GenTait(ρ, ρ0, c, refdensityfactor, exponent::Int)
    if exponent == 1
        return c^2 * (ρ - refdensityfactor * ρ0)
    else
        p0 = c^2 * ρ0 / exponent
        return p0 * ((ρ / ρ0)^exponent - refdensityfactor)
    end
end

"""
```julia
pressuretodensity_GenTait(p, ρ0, c, refdensityfactor, exponent)
```

Convert pressure to density

# Arguments:

  - `p::Float`:                particle pressure
  - `ρ0::Float`:               initial density
  - `c::Float`:                speed of sound
  - `refdensityfactor::Float`: reference density factor
  - `exponent::Int`:           exponent

# Returns:

  - particle density

# Example:

```jldoctest
# set equation of state
setequationsofstate(EquationOfState.GenTait)

# setup
c = 3.5
refdensityfactor = 0.9
p = 1.8
ρ0 = 0.4

# verify
exponent = 7
@assert pressuretodensity(p, ρ0, c, refdensityfactor, exponent) ≈ 0.477832244
exponent = 1
@assert pressuretodensity(p, ρ0, c, refdensityfactor, exponent) ≈ 0.506938775

# output

```
"""
function pressuretodensity_GenTait(p, ρ0, c, refdensityfactor, exponent::Int)
    if exponent == 1
        return p / c^2 + refdensityfactor * ρ0
    else
        p0 = c^2 * ρ0 / exponent
        return ρ0 * (p / p0 + refdensityfactor)^(1.0 / exponent)
    end
end

"""
```julia
densitytoenergy_GenTait(ρ, m, ρ0, c, refdensityfactor, exponent)
```

Convert density to energy

# Arguments:

  - `ρ::Float`:                particle density
  - `m::Float`:                particle mass
  - `ρ0::Float`:               initial density
  - `c::Float`:                speed of sound
  - `refdensityfactor::Float`: reference density factor
  - `exponent::Int`:           exponent

# Returns:

  - particle energy

# Example:

```jldoctest
# set equation of state
setequationsofstate(EquationOfState.GenTait)

# setup
c = 3.5
refdensityfactor = 0.9
m = 1.3
ρ0 = 0.4

# verify
exponent = 7
@assert densitytoenergy(ρ, m, ρ0, c, refdensityfactor, exponent) ≈ 0.1999413554e2
exponent = 1
@assert densitytoenergy(ρ, m, ρ0, c, refdensityfactor, exponent) ≈ -0.5553068732e1

# output

```
"""
function densitytoenergy_GenTait(ρ, m, ρ0, c, refdensityfactor, exponent::Int)
    if exponent == 1
        return -c^2 * m * (log(m^2 / (ρ0 * ρ)) - refdensityfactor * (1.0 + ρ0 / ρ))
    else
        p0 = c^2 * ρ0 / exponent
        return -p0 * (
            (1.0 / (1 - exponent)) * (m / (ρ0^exponent * ρ^(1 - exponent)) + m / ρ0) -
            refdensityfactor * (m / ρ0 + m / ρ)
        )
    end
end

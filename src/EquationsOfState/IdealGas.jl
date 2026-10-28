# ------------------------------------------------------------------------------
# Equations of state - ideal gas law
# ------------------------------------------------------------------------------

"""
```julia
densitytopressure_IdealGas(ρ, ρ0, c, refdensityfactor, exponent)
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
setequationsofstate(EquationOfState.IdealGas)

# setup
c = 3.5
ρ = 0.75

# verify
@assert densitytopressure(ρ, 0.0, c, 0.0, 0) ≈ 9.1875

# output

```
"""
function densitytopressure_IdealGas(ρ, ρ0, c, refdensityfactor, exponent::Int)
    return c^2 * ρ
end

"""
```julia
pressuretodensity_IdealGas(p, ρ0, c, refdensityfactor, exponent)
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
setequationsofstate(EquationOfState.IdealGas)

# setup
c = 3.5
p = 1.8
ρ0 = 0.4

# verify
@assert pressuretodensity(p, ρ0, c, 0.0, 0) ≈ 0.146938775

# output

```
"""
function pressuretodensity_IdealGas(p, ρ0, c, refdensityfactor, exponent::Int)
    return p / c^2
end

"""
```julia
densitytoenergy_IdealGas(ρ, m, ρ0, c, refdensityfactor, exponent)
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
setequationsofstate(EquationOfState.IdealGas)

# setup
c = 3.5
m = 1.3
ρ0 = 0.4

# verify
@assert densitytoenergy(ρ, m, ρ0, c, 0.0, 0) ≈ -0.2752956873e2

# output

```
"""
function densitytoenergy_IdealGas(ρ, m, ρ0, c, refdensityfactor, exponent::Int)
    return -c^2 * m * log(m^2 / (ρ0 * ρ))
end

# ------------------------------------------------------------------------------
# Momentum formulation - Adami
# ------------------------------------------------------------------------------

"""
```julia
specificcoefficient_Adami(ρ_i, ρ_j, m_i, m_j, dWdr_ij)
```

Compute the specific coefficient

# Arguments:

  - `ρ_i::Float`:     density of first particle
  - `ρ_j::Float`:     density of second particle
  - `m_i::Float`:     mass of first particle
  - `m_j::Float`:     mass of second particle
  - `dWdr_ij::Float`: first derivative of the kernel

# Returns:

  - specific coefficient

# Example:

```jldoctest
# setup
ρ_i = 1.1
ρ_j = 1.2
m_i = 2.2
m_j = 2.3
dWdr_ij = 3.3

# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# compute specific coefficient
speccoeff = specificcoefficient!(ρ_i, ρ_j, m_i, m_j, dWdr_ij)

# verify
@assert speccoeff ≈ 42.0

# output

```
"""
function specificcoefficient_Adami(ρ_i, ρ_j, m_i, m_j, dWdr_ij)
    return ((m_i / ρ_i)^2 + (m_j / ρ_j)^2) * (dWdr_ij / m_i)
end

"""
```julia
pressuregradient_Adami!(ρ_i, ρ_j, p_i, p_j, speccoeff, e_ij, a)
```

Compute the pressure gradient

# Arguments:

  - `ρ_i::Float`:       density of first particle
  - `ρ_j::Float`:       density of second particle
  - `p_i::Float`:       pressure of first particle
  - `p_j::Float`:       pressure of second particle
  - `speccoeff::Float`: specific coefficient
  - `e_ij::Array`:      normalized vector between particles
  - `a::Array`:         output acceleration array to sum in the pressure gradient

# Returns:

  - none

# Example:

```jldoctest
# setup
ρ_i = 1.1
ρ_j = 1.2
p_i = 2.2
p_j = 2.3
speccoeff = 3.3
e_ij = [4.0 / 5.0, 3.0 / 5.0]
a = [1.0, 2.0]

# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# compute pressure gradient
pressuregradient!(ρ_i, ρ_j, p_i, p_j, speccoeff, e_ij, a)

# verify
@assert acc ≈ [1.0, 2.0]

# output

```
"""
function pressuregradient_Adami!(ρ_i, ρ_j, p_i, p_j, speccoeff, e_ij, a)
    a += -(ρ_i * p_j + ρ_j * p_i) / (ρ_i + ρ_j) * speccoeff * e_ij
    return nothing
end

"""
```julia
shearforces_Adami!(
    ρ_i,
    ρ_j,
    v_i,
    v_j,
    kernelfac,
    μ_i,
    μ_j,
    μ_bulk_i,
    μ_bulk_j,
    absdist,
    speccoeff,
    e_ij,
    a,
)
```

Compute the shear forces

# Arguments:

  - `ρ_i::Float`:       density of first particle
  - `ρ_j::Float`:       density of second particle
  - `v_i::Array`:       velocity of first particle
  - `v_j::Array`:       velocity of second particle
  - `μ_i::Float`:       dynamic viscosity of first particle
  - `μ_j::Float`:       dynamic viscosity of second particle
  - `μ_bulk_i::Float`:  bulk viscosity of first particle
  - `μ_bulk_j::Float`:  bulk viscosity of second particle
  - `absdist::Float`:   absolute distance between particles
  - `speccoeff::Float`: specific coefficient
  - `e_ij::Array`:      normalized vector between particles
  - `a::Array`:         output acceleration array to sum in the shear forces

# Returns:

  - none

# Example:

```jldoctest
# setup
ρ_i = 1.1
ρ_j = 1.2
v_i = [1.0, 2.0]
v_j = [3.0, 4.0]
kernelfac = 2.2
μ_i = 3.3
μ_j = 3.4
μ_bulk_i = 4.4
μ_bulk_j = 4.5
absdist = 5.5
speccoef = 6.6
e_ij = [4.0 / 5.0, 3.0 / 5.0]
acc = [1.0, 2.0]

# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# compute shear forces
shearforces!(
    ρ_i,
    ρ_j,
    v_i,
    v_j,
    kernelfac,
    μ_i,
    μ_j,
    μ_bulk_i,
    μ_bulk_j,
    absdist,
    speccoeff,
    e_ij,
    a,
)

# verify
@assert acc ≈ [1.0, 2.0]

# output

```
"""
function shearforces_Adami!(
    ρ_i,
    ρ_j,
    v_i::Array,
    v_j::Array,
    kernelfac,
    μ_i,
    μ_j,
    μ_bulk_i,
    μ_bulk_j,
    absdist,
    speccoeff,
    e_ij::Array,
    a::Array,
)
    if !(μ_i > 0.0 && μ_j > 0.0)
        return
    end

    a += speccoeff * (2.0 * μ_i * μ_j / (μ_i + μ_j)) / absdist * (v_i - v_j)
    return nothing
end

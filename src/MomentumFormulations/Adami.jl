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
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# setup
ρ_i = 1.02
ρ_j = 0.97
m_i = 0.25
m_j = 0.23
dWdr_ij = 0.3

# compute specific coefficient
speccoeff = specificcoefficient(ρ_i, ρ_j, m_i, m_j, dWdr_ij)

# verify
@assert speccoeff ≈ ((m_i / ρ_i)^2 + (m_j / ρ_j)^2) * dWdr_ij / m_i

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
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# setup
ρ_i = 1.02
ρ_j = 0.97
p_i = 0.45
p_j = 0.43
speccoeff = 12.42
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
a = [1.0, 2.0, 3.0]

# compute pressure gradient
pressuregradient!(ρ_i, ρ_j, p_i, p_j, speccoeff, e_ij, a)

# verify
@assert a ≈ [1.0, 2.0, 3.0] - (ρ_i * p_j + ρ_j * p_i) / (ρ_i + ρ_j) * speccoeff * e_ij

# output

```
"""
function pressuregradient_Adami!(ρ_i, ρ_j, p_i, p_j, speccoeff, e_ij, a)
    a .+= -(ρ_i * p_j + ρ_j * p_i) / (ρ_i + ρ_j) * speccoeff * e_ij
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
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# setup
ρ_i = 1.02
ρ_j = 0.97
v_i = [0.1, 2.0, -0.4]
v_j = [0.27, -1.3, 0.5]
absdist = 0.3
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
speccoef = 12.42
kernelfac = 3
μ_i = 0.03
μ_j = 0.05
μ_bulk_i = 0.015
μ_bulk_j = 0.006
a = [1.0, 2.0, 3.0]

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
@assert a ≈
        [1.0, 2.0, 3.0] +
        speccoeff * (2.0 * μ_i * μ_j / (μ_i + μ_j)) / absdist * (v_i - v_j)

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

    a .+= speccoeff * (2.0 * μ_i * μ_j / (μ_i + μ_j)) / absdist * (v_i - v_j)
    return nothing
end

"""
```julia
standardbackgroundpressure_Adami!(ρ_i, ρ_j, p_background, speccoeff, e_ij, a_modified)
```

Compute the background pressure forces

# Arguments:

  - `ρ_i::Float`:          density of first particle
  - `ρ_j::Float`:          density of second particle
  - `p_background::Float`: background pressure
  - `speccoeff::Float`:    specific coefficient
  - `e_ij::Array`:         normalized vector between particles
  - `a_modified::Array`:   output modified acceleration array to sum in the background pressure forces

# Returns:

  - none

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# setup
ρ_i = 1.02
ρ_j = 0.97
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
speccoef = 12.42
p_background = 10.0
a_modified = [1.0, 2.0, 3.0]

# compute standard background pressure forces
standardbackgroundpressure!(ρ_i, ρ_j, p_background, speccoeff, e_ij, a_modified)

# verify
@assert a_modified ≈ [1.0, 2.0, 3.0] - speccoeff * p_background * e_ij

# output

```
"""
function standardbackgroundpressure_Adami!(
    ρ_i,
    ρ_j,
    p_background,
    speccoeff,
    e_ij::Array,
    a_modified::Array,
)
    a_modified .+= -speccoeff * p_background * e_ij
    return nothing
end

"""
```julia
generalizedbackgroundpressure_Adami!(
    ρ_i,
    ρ_j,
    m_i,
    m_j,
    p_background_mod,
    speccoeff,
    e_ij,
    a_modified,
)
```

Compute the generalized background pressure forces

# Arguments:

  - `ρ_i::Float`:              density of first particle
  - `ρ_j::Float`:              density of second particle
  - `m_i::Float`:              mass of first particle
  - `m_j::Float`:              mass of second particle
  - `p_background_mod::Float`: modified background pressure
  - `dWdr_ij_mod::Float`:      modified derivative of kernel factor
  - `e_ij::Array`:             normalized vector between particles
  - `a_modified::Array`:       output modified acceleration array to sum in the background pressure forces

# Returns:

  - none

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# setup
ρ_i = 1.02
ρ_j = 0.97
m_i = 0.25
m_j = 0.23
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
p_background_mod = 8.4
dWdr_ij_mod = 0.32
a_modified = [1.0, 2.0, 3.0]

# compute generalized background pressure forces
generalizedbackgroundpressure!(
    ρ_i,
    ρ_j,
    m_i,
    m_j,
    p_background_mod,
    dWdr_ij_mod,
    e_ij,
    a_modified,
)

# verify
@assert a_modified ≈ [1.0, 2.0, 3.0] - p_background_mod * m_i / ρ_i^2 * dWdr_ij_mod * e_ij

# output

```
"""
function generalizedbackgroundpressure_Adami!(
    ρ_i,
    ρ_j,
    m_i,
    m_j,
    p_background_mod,
    dWdr_ij_mod,
    e_ij::Array,
    a_modified::Array,
)
    a_modified .+= -p_background_mod * m_i / ρ_i^2 * dWdr_ij_mod * e_ij
    return nothing
end

"""
```julia
modifiedvelocitycontribution_Adami!(
    ρ_i,
    ρ_j,
    v_i,
    v_j,
    v_modified_i,
    v_modified_j,
    speccoeff,
    e_ij,
    a_modified,
)
```

Compute the modified velocity force contribution

# Arguments:

  - `ρ_i::Float`:          density of first particle
  - `ρ_j::Float`:          density of second particle
  - `v_i::Array`:          velocity of first particle
  - `v_j::Array`:          velocity of second particle
  - `v_modified_i::Array`: modified velocity of first particle
  - `v_modified_j::Array`: modified velocity of second particle
  - `speccoeff::Float`:    specific coefficient
  - `e_ij::Array`:         normalized vector between particles
  - `a_modified::Array`:   output modified acceleration array to sum in the background pressure forces

# Returns:

  - none

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# setup
ρ_i = 1.02
ρ_j = 0.97
v_i = [0.1, 2.0, -0.4]
v_j = [0.27, -1.3, 0.5]
v_modified_i = [0.12, 1.9, -0.42]
v_modified_j = [0.3, -1.2, 0.77]
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
speccoeff = 12.42
a_modified = [1.0, 2.0, 3.0]

# compute modified velocity force contribution
modifiedvelocitycontribution!(
    ρ_i,
    ρ_j,
    v_i,
    v_j,
    v_modified_i,
    v_modified_j,
    speccoeff,
    e_ij,
    a_modified,
)

# verify
A_i = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
A_j = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
for i = 1:3
    for j = 1:3
        A_i[i, j] = ρ_i * v_i[i] * (v_modified_i[j] - v_i[j])
        A_j[i, j] = ρ_j * v_j[i] * (v_modified_j[j] - v_j[j])
    end
end
A_ij_e_ij = 0.5 * A_i * e_ij + 0.5 * A_j * e_ij

@assert a_modified ≈ [1.0, 2.0, 3.0] + speccoef * A_ij_e_ij

# output

```
"""
function modifiedvelocitycontribution_Adami!(
    ρ_i,
    ρ_j,
    v_i::Array,
    v_j::Array,
    v_modified_i::Array,
    v_modified_j::Array,
    speccoeff,
    e_ij::Array,
    a_modified::Array,
)
    A_ij_e_ij = [0.0, 0.0, 0.0]

    if true
        A_ij_e_ij .+= 0.5 * ρ_i * ((v_modified_i - v_i) ⋅ e_ij) * v_i
    end
    if true
        A_ij_e_ij .+= 0.5 * ρ_j * ((v_modified_j - v_j) ⋅ e_ij) * v_j
    end

    a_modified .+= speccoeff * A_ij_e_ij
    return nothing
end

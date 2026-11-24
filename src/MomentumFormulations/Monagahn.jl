"""
```julia
specificcoefficient_Monagahn(ρ_i, ρ_j, m_i, m_j, dWdr_ij)
```

Compute the specific coefficient

# Arguments:

  - `ρ_i::PulvisScalar`:     density of first particle
  - `ρ_j::PulvisScalar`:     density of second particle
  - `m_i::PulvisScalar`:     mass of first particle
  - `m_j::PulvisScalar`:     mass of second particle
  - `dWdr_ij::PulvisScalar`: first derivative of the kernel

# Returns:

  - specific coefficient

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Monagahn)

# setup
ρ_i = 1.02
ρ_j = 0.97
m_i = 0.25
m_j = 0.23
dWdr_ij = 0.3

# compute specific coefficient
speccoeff = Pulvis.specificcoefficient(ρ_i, ρ_j, m_i, m_j, dWdr_ij)

# verify
@assert speccoeff ≈ dWdr_ij * m_i

# output

```
"""
function specificcoefficient_Monagahn(
    ρ_i::PulvisScalar,
    ρ_j::PulvisScalar,
    m_i::PulvisScalar,
    m_j::PulvisScalar,
    dWdr_ij::PulvisScalar,
)
    return dWdr_ij * m_i
end

"""
```julia
pressuregradient_Monagahn!(ρ_i, ρ_j, p_i, p_j, speccoeff, e_ij, a)
```

Compute the pressure gradient

# Arguments:

  - `ρ_i::PulvisScalar`:         density of first particle
  - `ρ_j::PulvisScalar`:         density of second particle
  - `p_i::PulvisScalar`:         pressure of first particle
  - `p_j::PulvisScalar`:         pressure of second particle
  - `speccoeff::PulvisScalar`:   specific coefficient
  - `e_ij::Array{PulvisScalar}`: normalized vector between particles
  - `a::Array{PulvisScalar}`:    output acceleration array to sum in the pressure gradient

# Returns:

  - nothing

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Monagahn)

# setup
ρ_i = 1.02
ρ_j = 0.97
p_i = 0.45
p_j = 0.43
speccoeff = 12.24
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
a = [1.0, 2.0, 3.0]

# compute pressure gradient
Pulvis.pressuregradient!(ρ_i, ρ_j, p_i, p_j, speccoeff, e_ij, a)

# verify
@assert a ≈ [1.0, 2.0, 3.0] - (p_i / ρ_i^2 + p_j / ρ_j^2) * speccoeff * e_ij

# output

```
"""
function pressuregradient_Monagahn!(
    ρ_i::PulvisScalar,
    ρ_j::PulvisScalar,
    p_i::PulvisScalar,
    p_j::PulvisScalar,
    speccoeff::PulvisScalar,
    e_ij::AbstractArray{PulvisScalar},
    a::AbstractArray{PulvisScalar},
)
    a .+= -(p_i / ρ_i^2 + p_j / ρ_j^2) * speccoeff * e_ij
    return nothing
end

"""
```julia
shearforces_Monagahn!(
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

  - `ρ_i::PulvisScalar`:         density of first particle
  - `ρ_j::PulvisScalar`:         density of second particle
  - `v_i::Array{PulvisScalar}`:  velocity of first particle
  - `v_j::Array{PulvisScalar}`:  velocity of second particle
  - `kernelfac::Int`:            kernel factor
  - `μ_i::PulvisScalar`:         dynamic viscosity of first particle
  - `μ_j::PulvisScalar`:         dynamic viscosity of second particle
  - `μ_bulk_i::PulvisScalar`:    bulk viscosity of first particle
  - `μ_bulk_j::PulvisScalar`:    bulk viscosity of second particle
  - `absdist::PulvisScalar`:     absolute distance between particles
  - `speccoeff::PulvisScalar`:   specific coefficient
  - `e_ij::Array{PulvisScalar}`: normalized vector between particles
  - `a::Array{PulvisScalar}`:    output acceleration array to sum in the shear forces

# Returns:

  - nothing

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Monagahn)

# setup
ρ_i = 1.02
ρ_j = 0.97
v_i = [0.1, 2.0, -0.4]
v_j = [0.27, -1.3, 0.5]
absdist = 0.3
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
speccoeff = 12.42
kernelfac = 3
μ_i = 0.03
μ_j = 0.05
μ_bulk_i = 0.015
μ_bulk_j = 0.006
a = [1.0, 2.0, 3.0]

# compute shear forces
Pulvis.shearforces!(
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
μ_scaled = (μ_i > 0.0 && μ_j > 0.0) ? (2.0 * μ_i * μ_j / (3.0 * (μ_i + μ_j))) : 0.0
μ_bulk =
    (μ_bulk_i > 0.0 && μ_bulk_j > 0.0) ?
    (2.0 * μ_bulk_i * μ_bulk_j / (μ_bulk_i + μ_bulk_j)) : 0.0
v_ij = (v_i - v_j)
inv_ρ_factor = 1.0 / (ρ_i * ρ_j * absdist)
diffusion_coeff = 5.0 * μ_scaled - μ_bulk
convection_coeff = kernelfac * (μ_bulk + μ_scaled)
@assert a ≈
        [1.0, 2.0, 3.0] +
        speccoeff * (
    diffusion_coeff * inv_ρ_factor * v_ij +
    convection_coeff * (v_ij ⋅ e_ij) * inv_ρ_factor * e_ij
)

# output

```
"""
function shearforces_Monagahn!(
    ρ_i::PulvisScalar,
    ρ_j::PulvisScalar,
    v_i::AbstractArray{PulvisScalar},
    v_j::AbstractArray{PulvisScalar},
    kernelfac::Int,
    μ_i::PulvisScalar,
    μ_j::PulvisScalar,
    μ_bulk_i::PulvisScalar,
    μ_bulk_j::PulvisScalar,
    absdist::PulvisScalar,
    speccoeff::PulvisScalar,
    e_ij::AbstractArray{PulvisScalar},
    a::AbstractArray{PulvisScalar},
)
    μ_scaled = (μ_i > 0.0 && μ_j > 0.0) ? (2.0 * μ_i * μ_j / (3.0 * (μ_i + μ_j))) : 0.0
    μ_bulk =
        (μ_bulk_i > 0.0 && μ_bulk_j > 0.0) ?
        (2.0 * μ_bulk_i * μ_bulk_j / (μ_bulk_i + μ_bulk_j)) : 0.0
    inv_ρ_factor = 1.0 / (ρ_i * ρ_j * absdist)
    v_ij = (v_i - v_j)

    # diffusion
    diffusion_coeff = 5.0 * μ_scaled - μ_bulk
    a .+= speccoeff * diffusion_coeff * inv_ρ_factor * v_ij

    # convection
    convection_coeff = kernelfac * (μ_bulk + μ_scaled)
    a .+= speccoeff * convection_coeff * (v_ij ⋅ e_ij) * inv_ρ_factor * e_ij
    return nothing
end

"""
```julia
standardbackgroundpressure_Monagahn!(ρ_i, ρ_j, p_background, speccoeff, e_ij, a_modified)
```

Compute the background pressure forces

# Arguments:

  - `ρ_i::PulvisScalar`:               density of first particle
  - `ρ_j::PulvisScalar`:               density of second particle
  - `p_background::PulvisScalar`:      background pressure
  - `speccoeff::PulvisScalar`:         specific coefficient
  - `e_ij::Array{PulvisScalar}`:       normalized vector between particles
  - `a_modified::Array{PulvisScalar}`: output modified acceleration array to sum in the background pressure forces

# Returns:

  - nothing

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Monagahn)

# setup
ρ_i = 1.02
ρ_j = 0.97
e_ij = [1.0 / √21.0, 2.0 / √21.0, 4.0 / √21.0]
speccoeff = 12.42
p_background = 10.0
a_modified = [1.0, 2.0, 3.0]

# compute standard background pressure forces
Pulvis.standardbackgroundpressure!(ρ_i, ρ_j, p_background, speccoeff, e_ij, a_modified)

# verify
@assert a_modified ≈
        [1.0, 2.0, 3.0] - speccoeff * p_background * (1.0 / ρ_i^2 + 1.0 / ρ_j^2) * e_ij

# output

```
"""
function standardbackgroundpressure_Monagahn!(
    ρ_i::PulvisScalar,
    ρ_j::PulvisScalar,
    p_background::PulvisScalar,
    speccoeff::PulvisScalar,
    e_ij::AbstractArray{PulvisScalar},
    a_modified::AbstractArray{PulvisScalar},
)
    a_modified .+= -speccoeff * p_background * (1.0 / ρ_i^2 + 1.0 / ρ_j^2) * e_ij
    return nothing
end

"""
```julia
generalizedbackgroundpressure_Monagahn!(
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

  - `ρ_i::PulvisScalar`:               density of first particle
  - `ρ_j::PulvisScalar`:               density of second particle
  - `m_i::PulvisScalar`:               mass of first particle
  - `m_j::PulvisScalar`:               mass of second particle
  - `p_background_mod::PulvisScalar`:  modified background pressure
  - `dWdr_ij_mod::PulvisScalar`:       modified derivative of kernel factor
  - `e_ij::Array{PulvisScalar}`:       normalized vector between particles
  - `a_modified::Array{PulvisScalar}`: output modified acceleration array to sum in the background pressure forces

# Returns:

  - nothing

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Monagahn)

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
Pulvis.generalizedbackgroundpressure!(
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
@assert a_modified ≈ [1.0, 2.0, 3.0] - p_background_mod * (m_j / ρ_i^2) * dWdr_ij_mod * e_ij

# output

```
"""
function generalizedbackgroundpressure_Monagahn!(
    ρ_i::PulvisScalar,
    ρ_j::PulvisScalar,
    m_i::PulvisScalar,
    m_j::PulvisScalar,
    p_background_mod::PulvisScalar,
    dWdr_ij_mod::PulvisScalar,
    e_ij::AbstractArray{PulvisScalar},
    a_modified::AbstractArray{PulvisScalar},
)
    a_modified .+= -p_background_mod * m_j / ρ_i^2 * dWdr_ij_mod * e_ij
    return nothing
end

"""
```julia
modifiedvelocitycontribution_Monagahn!(
    ρ_i,
    ρ_j,
    v_i,
    v_j,
    speccoeff,
    e_ij,
    a_modified;
    v_modified_i,
    v_modified_j,
)
```

Compute the modified velocity force contribution

# Arguments:

  - `ρ_i::PulvisScalar`:                 density of first particle
  - `ρ_j::PulvisScalar`:                 density of second particle
  - `v_i::Array{PulvisScalar}`:          velocity of first particle
  - `v_j::Array{PulvisScalar}`:          velocity of second particle
  - `speccoeff::PulvisScalar`:           specific coefficient
  - `e_ij::Array{PulvisScalar}`:         normalized vector between particles
  - `a_modified::Array{PulvisScalar}`:   output modified acceleration array to sum in the background pressure forces
  - `v_modified_i::Array{PulvisScalar}`: modified velocity of first particle, optional
  - `v_modified_j::Array{PulvisScalar}`: modified velocity of second particle, optional

# Returns:

  - nothing

# Example:

```jldoctest
# set the momentum formulation
setmomentumformulation(MomentumFormulation.Monagahn)

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
Pulvis.modifiedvelocitycontribution!(
    ρ_i,
    ρ_j,
    v_i,
    v_j,
    speccoeff,
    e_ij,
    a_modified,
    v_modified_i = v_modified_i,
    v_modified_j = v_modified_j,
)

# verify
A_i = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
A_j = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
for i ∈ 1:3
    for j ∈ 1:3
        A_i[i, j] = ρ_i * v_i[i] * (v_modified_i[j] - v_i[j])
        A_j[i, j] = ρ_j * v_j[i] * (v_modified_j[j] - v_j[j])
    end
end
A_ij_e_ij = 1.0 / ρ_i^2 * A_i * e_ij + 1.0 / ρ_j^2 * A_j * e_ij

@assert a_modified ≈ [1.0, 2.0, 3.0] + speccoeff * A_ij_e_ij

# try with no v_modified_*
a_modified = [1.0, 2.0, 3.0]
Pulvis.modifiedvelocitycontribution!(ρ_i, ρ_j, v_i, v_j, speccoeff, e_ij, a_modified)

# verify
@assert a_modified ≈ [1.0, 2.0, 3.0]

# output

```
"""
function modifiedvelocitycontribution_Monagahn!(
    ρ_i::PulvisScalar,
    ρ_j::PulvisScalar,
    v_i::AbstractArray{PulvisScalar},
    v_j::AbstractArray{PulvisScalar},
    speccoeff::PulvisScalar,
    e_ij::AbstractArray{PulvisScalar},
    a_modified::AbstractArray{PulvisScalar};
    v_modified_i = nothing,
    v_modified_j = nothing,
)
    A_ij_e_ij = [0.0, 0.0, 0.0]

    if v_modified_i != nothing
        A_ij_e_ij .+= ((v_modified_i - v_i) ⋅ e_ij) / ρ_i * v_i
    end
    if v_modified_j != nothing
        A_ij_e_ij .+= ((v_modified_j - v_j) ⋅ e_ij) / ρ_j * v_j
    end

    a_modified .+= speccoeff * A_ij_e_ij
    return nothing
end

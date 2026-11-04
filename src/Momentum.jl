# ------------------------------------------------------------------------------
# Momentum computation at a single point
# ------------------------------------------------------------------------------

"""
```julia
momentumequationneighbors!(
    particleindex,
    particleids,
    typeids,
    neighborindexes,
    μ_dynamic,
    μ_bulk,
    dimension,
    r,
    radius,
    m,
    ρ,
    p,
    v,
    a,
    v_modified,
    a_modified,
)
```

Sum the momentum contributions to the particle from neighbor particles

# Arguments:

  - `particleindex::Int`:     global index number of current particle
  - `particleids::Array`:     input array of particle ids by type
  - `neighborindexes::Array`: input array of neighbor particle global indexes
  - `typeids::Array`:         input array of particle type ids
  - `μ_dynamic::Array`        input array of dynamic viscosity for each particle type
  - `μ_bulk::Array`           input array of bulk viscosity for each particle type
  - `dimension::Int`:         spatial dimension of problem
  - `r::StateData`:           input arrays of particle positions by particle type and id
  - `radius::StateData`:      input arrays of particle radi by particle type and id
  - `m::StateData`:           input arrays of particle masses by particle type and id
  - `ρ::StateData`:           input arrays of particle densities by particle type and id
  - `p::StateData`:           input arrays of particle pressures by particle type and id
  - `v::StateData`:           input arrays of particle velocities by particle type and id
  - `a::StateData`:           output arrays of particle acceleration by particle type and id
  - `v_modified::StateData`:  output arrays of modified particle velocities by particle type and id
  - `a_modified::StateData`:  output arrays of modified particle acceleration by particle type and id

# Returns:

  - none

# Example:

```jldoctest
# setup
particleindex = 1
particleids = [1, 2, 1]
typeids = [1, 1, 2]
neighborindexes = [3, 5, 2, 3]
μ_dynamic = [1.1, 1.2]
μ_bulk = [2.1, 2.2]
dimension = 2
r = StateData{2}([[1.0, 1.0, 1.1, 1.2], [1.3, 1.4]])
radius = StateData{2}([[1.0, 1.1], [1.2]])
m = StateData{2}([[42.0, 13.0], [14.0]])
ρ = StateData{2}([[1.0, 2.0], [3.0]])
p = StateData{2}([[1.0, 2.0], [3.0]])
v = StateData{2}([[1.0, 2.0, 1.1, 2.1], [3.0, 3.1]])
a = StateData{2}([[1.0, 2.0, 0.0, 0.0], [0.0, 0.0]])
v_modified = StateData{2}([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0]])
a_modified = StateData{2}([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0]])

# set shape function
setkernelshape(Kernel.CubicSpline)

# set the momentum formulation
setmomentumformulation(MomentumFormulation.Monagahn)

# compute density
momentumequationneighbors!(
    particleindex,
    particleids,
    typeids,
    neighborindexes,
    μ_dynamic,
    μ_bulk,
    dimension,
    r,
    radius,
    m,
    ρ,
    p,
    v,
    a,
    v_modified,
    a_modified,
)

# verify
@assert a[typeids[particleindex]][((particleids[particleindex]-1)*dimension+1):(particleids[particleindex]*dimension)] ≈
        [2778.2561852752383, 3851.6790419419376]
neighborindex = neighborindexes[3]
@assert a[typeids[neighborindex]][((particleids[neighborindex]-1)*dimension+1):(particleids[neighborindex]*dimension)] ≈
        [0.0, 0.0]
neighborindex = neighborindexes[4]
@assert a[typeids[neighborindex]][((particleids[neighborindex]-1)*dimension+1):(particleids[neighborindex]*dimension)] ≈
        [0.0, 0.0]

# output

```
"""
function momentumequationneighbors!(
    particleindex::Int,
    particleids::Array,
    typeids::Array,
    neighborindexes::Array,
    μ_dynamic::Array,
    μ_bulk::Array,
    dimension::Int,
    r::StateData,
    radius::StateData,
    m::StateData,
    ρ::StateData,
    p::StateData,
    v::StateData,
    a::StateData,
    v_modified::StateData,
    a_modified::StateData,
)
    particleid = particleids[particleindex]
    typeid = typeids[particleindex]

    for i = neighborindexes[particleindex]:(neighborindexes[particleindex+1]-1)
        neighborindex = neighborindexes[i]
        neighborid = particleids[neighborindex]
        neighbortypeid = typeids[neighborindex]

        r_ij =
            view(r[neighbortypeid], ((neighborid-1)*dimension+1):(neighborid*dimension)) -
            view(r[typeid], ((particleid-1)*dimension+1):(particleid*dimension))
        absdistance = norm(r_ij)

        # specific coefficient
        dWdr_ij = dwdrij(absdistance, radius[typeid][particleid], dimension)

        speccoeff = specificcoefficient(
            ρ[typeid][particleid],
            ρ[neighbortypeid][neighborid],
            m[typeid][particleid],
            m[neighbortypeid][neighborid],
            dWdr_ij,
        )

        # pressure gradient
        e_ij = (-1.0 / absdistance) * r_ij

        pressuregradient!(
            ρ[typeid][particleid],
            ρ[neighbortypeid][neighborid],
            p[typeid][particleid],
            p[neighbortypeid][neighborid],
            speccoeff,
            e_ij,
            view(a[typeid], ((particleid-1)*dimension+1):(particleid*dimension)),
        )

        # shear forces
        shearforces!(
            ρ[typeid][particleid],
            ρ[typeid][neighborid],
            view(v[typeid], ((particleid-1)*dimension+1):(particleid*dimension)),
            view(v[neighbortypeid], ((neighborid-1)*dimension+1):(neighborid*dimension)),
            dimension + 2,
            μ_dynamic[typeid],
            μ_dynamic[neighbortypeid],
            μ_bulk[typeid],
            μ_bulk[neighbortypeid],
            absdistance,
            speccoeff,
            e_ij,
            view(a[typeid], ((particleid-1)*dimension+1):(particleid*dimension)),
        )

        # transport velocity - TODO

        # artificial viscosity - TODO
    end
    return nothing
end

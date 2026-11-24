# ------------------------------------------------------------------------------
# Smoothed Particle Hydrodynamics
# ------------------------------------------------------------------------------

"""
```julia
preevaluatetimestep!(
    bininfo,
    dimension,
    time,
    typeinfo,
    stateinfo,
    statedata,
    particleids,
    typeids,
    neighborindexes,
)
```

Prepare for timestep; apply the open boundary conditions, update local particle id and particle type id arrays and the neighbor particle global index array

# Arguments:

  - `bininfo::BinInfo`:            info about background bin layout
  - `time::PulvisScalar`:          current solver time
  - `dimension::Int`:              spatial dimension of problem
  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays
  - `particleids::Array{Int}`:     input/output array of local particle ids by type
  - `typeids::Array{Int}`:         input/output array of particle type ids
  - `neighborindexes::Array{Int}`: input/output array of neighbor particle global indexes

# Returns:

  - nothing
"""
function preevaluatetimestep_SPH!(
    bininfo::BinInfo,
    dimension::Int,
    time::PulvisScalar,
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
    particleids::Array{Int},
    typeids::Array{Int},
    neighborindexes::Array{Int},
)
    # Get state data
    boundaryid = getparticlestatedata(statedata, "boundary id")
    r = getparticlestatedata(statedata, "position")
    v = getparticlestatedata(statedata, "velocity")
    radius = getparticlestatedata(statedata, "radius")

    # Index of open boundary type
    boundaryindex = findfirst(t -> t.name == "open boundary", typeinfo)

    # Boundaries
    if boundaryindex != nothing
        prescribeopenboundarystates!(
            time,
            dimension,
            boundaryid[boundaryindex],
            r[boundaryindex],
            v[boundaryindex],
        )
    end

    # Build pair list
    buildpairlist!(
        typeinfo,
        bininfo,
        dimension,
        r,
        radius,
        particleids,
        typeids,
        nighborindexes,
    )
    return nothing
end

"""
```julia
evaluatetimestep!(
    particleids,
    typeids,
    neighborindexes,
    dimension,
    typeinfo,
    stateinfo,
    statedata,
)
```

Evaluate time step; compute density, convert density to pressure, then various forces, acceleration, etc from density and pressure

# Arguments:

  - `particleids::Array{Int}`:     input array of local particle ids by type
  - `typeids::Array{Int}`:         input array of particle type ids
  - `neighborindexes::Array{Int}`: input array of neighbor particle global indexes
  - `dimension::Int`:              spatial dimension of problem
  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing
"""
function evaluatetimestep_SPH!(
    particleids::Array{Int},
    typeids::Array{Int},
    neighborindexes::Array{Int},
    dimension::Int,
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    # Get state data
    r = getparticlestatedata(statedata, "position")
    radius = getparticlestatedata(statedata, "radius")
    m = getparticlestatedata(statedata, "mass")
    ρ = getparticlestatedata(statedata, "density")
    p = getparticlestatedata(statedata, "pressure")
    v = getparticlestatedata(statedata, "velocity")
    a = getparticlestatedata(statedata, "acceleration")
    v_modified = getparticlestatedata(statedata, "modified velocity")
    a_modified = getparticlestatedata(statedata, "modified acceleration")

    # Compute density field
    computedensity!(particleids, typeids, neighborindexes, dimension, r, radius, m, ρ)

    # Convert from density to pressure
    computepressure!(particleids, typeids, ρ, ρ0, c, refdensityfactor, exponent, p)

    # Compute various forces, acceleration, etc from density and pressure
    computeforces!(
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
    return nothing
end

"""
```julia
postevaluatetimestep!(
    particleids,
    typeids,
    neighborindexes,
    dimension,
    typeinfo,
    stateinfo,
    statedata,
)
```

Conclude time step; check boundary phase changes, evaluate peridynamic damage, etc

# Arguments:

  - `maxinteractiondistance::PulvisScalar`: max distance between particles for interaction
  - `dimension::Int`:                       spatial dimension of problem
  - `typeinfo::Array{TypeInfo}`:            input/output array of particle type information
  - `stateinfo::Array{StateInfo}`:          input array of particle state information
  - `statedata::Dict`:                      input/output dictionary of particle state data arrays

# Returns:

  - nothing
"""
function postevaluatetimestep_SPH!(
    maxinteractiondistance::PulvisScalar,
    dimension::Int,
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    # Boundaries
    checkopenboundaryphasechange!(
        maxinteractiondistance,
        dimension,
        typeinfo,
        stateinfo,
        statedata,
    )

    # Standard phase change - TODO
    # evaluatephasechange!(typeinfo, stateinfo, statedata)

    # Peridynamics - TODO
    # evaluatedamage(typeinfo, stateinfo, statedata)
    return nothing
end

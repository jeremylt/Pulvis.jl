# ------------------------------------------------------------------------------
# Open boundary conditions
# ------------------------------------------------------------------------------

boundaryfunctions = []
boundarynormal::Array{Array{PulvisScalar}} = [[]]

"""
```julia
prescribeopenboundarystates!(time, dimension, boundaryid, r, v)
```

Set values of boundary particles based upon prescribed boundary functions

# Arguments:

  - `time::PulvisScalar`:     current solver time
  - `dimension::Int`:         spatial dimension of problem
  - `boundaryid::Array{Int}`: input array of particle boundary ids
  - `r::Array{PulvisScalar}`: input array of particle positions
  - `v::Array{PulvisScalar}`: output array of particle velocities

# Returns:

  - none
"""
function prescribeopenboundarystates!(
    time::PulvisScalar,
    dimension::Int,
    boundaryid::Array{Int},
    r::Array{PulvisScalar},
    v::Array{PulvisScalar},
)
    numparticles = length(boundaryid)

    for particleid_i ∈ 1:numparticles
        boundaryid_i = boundaryid[particleid_i]
        r_i = view(r, ((particleid_i-1)*dimension+1):(particleid_i*dimension))

        v[((particleid_i-1)*dimension+1):(particleid_i*dimension)] .=
            -boundaryfunctions[boundaryid_i](r_i, time, 0) * boundarynormal[boundaryid_i]
    end
    return nothing
end

"""
```julia
checkopenboundaryphasechange!(maxinteractiondistance, dimension, typeinfo, stateinfo, statedata)
```

Set values of boundary particles based upon prescribed boundary functions

# Arguments:

  - `maxinteractiondistance::PulvisScalar`: max distance between particles for interaction
  - `dimension::Int`:                       spatial dimension of problem
  - `typeinfo::Array{TypeInfo}`:            input/output array of particle type information
  - `stateinfo::Array{StateInfo}`:          input array of particle state information
  - `statedata::Dict`:                      input/output dictionary of particle state data arrays

# Returns:

  - none
"""
function checkopenboundaryphasechange!(
    maxinteractiondistance::PulvisScalar,
    dimension::Int,
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    # Parameters
    numparticlesperdirection = Int(round(maxinteractiondistance ./ initialparticlespacing))
    toleranceopenboundarytofluid = 0.1 * initialparticlespacing

    # Particles to move and remove
    particlestomovetofluid::Vector{Int} = []
    particlestoremoveboundary::Vector{Int} = []
    particlestomovetoboundary::Vector{Int} = []
    particlestoremovefluid::Vector{Int} = []

    # Boundary particles
    boundarytypeindex = findfirst(t -> t.name == "open boundary", typeinfo)
    if boundarytypeindex != nothing
        boundaryid = getparticlestatedata(statedata, "boundary id")[boundarytypeindex]
        r = getparticlestatedata(statedata, "position")[boundarytypeindex]
        numparticles = typeinfo[boundarytypeindex].numparticles

        for particleid_i ∈ 1:numparticles
            boundaryid_i = boundaryid[particleid_i]
            r_i = view(r, ((particleid_i-1)*dimension+1):(particleid_i*dimension))

            distancefromplane =
                boundarynormal[boundaryid_i] ⋅ (r_i .- planepoints[boundaryid_i])
            if distancefromplane < -toleranceopenboundarytofluid[boundaryid_i]
                # Change over to fluid phase
                push!(particlestomovetofluid, particleid_i)
                r_i .+=
                    numparticlesperdirection[boundaryid_i] *
                    boundaryinitialparticlespacing[boundaryid_id] *
                    boundarynormal[boundaryid_i]
            elseif distancefromplane > maxinteractiondistance
                push!(particlestoremove, particleid_i)
            end
        end
    end

    # Fluid particles
    fluidtypeindex = findfirst(t -> t.name == "fluid phase", typeinfo)
    if fluidtypeindex != nothing
        r = getparticlestatedata(statedata, "position")[fluidtypeindex]
        numparticles = typeinfo[fluidtypeindex].numparticles
        numboundaries = length(boundarynormals)

        for particleid_i ∈ 1:numparticles
            for boundaryid_j ∈ 1:numboundaries
                distancefromplane =
                    boundarynormal[boundaryid_j] ⋅ (r_i .- planepoints[boundaryid_j])
                if distancefromplane > 0.0
                    push!(particlestomovetoboundary, particleid_i)
                end
            end
        end
    end

    # Move boundary particles
    copyparticles!(
        "open boundary",
        "fluid phase",
        particlestomovetofluid,
        particlestoremoveboundary,
        typeinfo,
        stateinfo,
        statedata,
    )
    removeparticles!(
        "open boundary",
        particlestoremoveboundary,
        typeinfo,
        stateinfo,
        statedata,
    )
    # Move fluid particles
    copyparticles!(
        "fluid phase",
        "open boundary",
        particlestomovetoboundary,
        particlestoremovefluid,
        typeinfo,
        stateinfo,
        statedata,
    )
    removeparticles!("fluid phase", particlestoremovefluid, typeinfo, stateinfo, statedata)
end

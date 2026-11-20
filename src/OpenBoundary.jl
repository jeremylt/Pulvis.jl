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
    numparticles = Int(length(r) / dimension)

    for particleid ∈ 1:numparticles
        boundaryid_i = boundaryid[particleid]
        r_i = view(r, ((particleid-1)*dimension+1):(particleid*dimension))

        v[((particleid-1)*dimension+1):(particleid*dimension)] .=
            -boundaryfunctions[boundaryid_i](r_i, time, 0) * boundarynormal[boundaryid]
    end
    return nothing
end

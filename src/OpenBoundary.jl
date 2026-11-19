# ------------------------------------------------------------------------------
# Open boundary conditions
# ------------------------------------------------------------------------------

function prescribeopenboundarystates!(
    time::PulvisScalar,
    dimension::Int,
    boundaryid::Array{Int},
    r::Array{PulvisScalar},
    v::Array{PulvisScalar},
)
    numparticles = Int(length(r) / dimension)

    for particleid = 1:numparticles
        boundaryid_i = boundaryid[particleid]
        r_i = view(r, ((particleid-1)*dimension+1):(particleid*dimension))

        v[((particleid-1)*dimension+1):(particleid*dimension)] .=
            -boundaryfunction![boundaryid_i](r_i, time, 0) * boundarynormal[boundrayid]
    end
    return nothing
end

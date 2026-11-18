# ------------------------------------------------------------------------------
# Time step hooks
# ------------------------------------------------------------------------------

function pretimestep(
    particleids::Array{Int},
    typeids::Array{Int},
    neighborindexes::Array{Int},
    bininfo::BinInfo,
    dimension::Int,
    time::PulvisScalar,
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    # Boundaries
    perscribeopenboundarystates!(time, typeinfo, stateinfo, statedata)

    # Get state data
    r = getparticlestatedata(statedata, "position")
    radius = getparticlestatedata(statedata, "radius")

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
end

function evaluatetimestep(
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
end

function posttimestep(typeinfo::Array, stateinfo::Array, statedata::Dict)
    # Boundaries
    checkopenboundaryphasechange(time, typeinfo, stateinfo, statedata)
    evaluatephasechange(time, typeinfo, stateinfo, statedata)

    # Peridynamics
    # evaluatedamage(time, typeinfo, stateinfo, statedata)
end

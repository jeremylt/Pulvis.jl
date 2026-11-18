# ------------------------------------------------------------------------------
# Time step hooks
# ------------------------------------------------------------------------------

function pretimestep(time, typeinfo::Array, stateinfo::Array, statedata::Dict)
    # Boundaries
    perscribeopenboundarystates(time, typeinfo, stateinfo, statedata)
end

function evaluatetimestep(
    dimension::Int,
    typeinfo::Array,
    stateinfo::Array,
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

    # Build pair list
    buildpairlist!(particleids, typeids, dimension, r, radius, nighborindexes)

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

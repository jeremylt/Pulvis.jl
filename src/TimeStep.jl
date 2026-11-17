# ------------------------------------------------------------------------------
# Time step hooks
# ------------------------------------------------------------------------------

function pretimestep(time, typeinfo::Array, stateinfo::Array, statedata::Dict)
    # Boundaries
    perscribeopenboundarystates(time, typeinfo, stateinfo, statedata)
end

function evaluatetimestep(typeinfo::Array, stateinfo::Array, statedata::Dict)
    # Setup for time step
    buildpairlist!()

    # Compute density field
    computedensity!()

    # Convert from density to pressure
    computepressure!()

    # Compute various forces, acceleration, etc from density and pressure
    computeforces!()
end

function posttimestep(typeinfo::Array, stateinfo::Array, statedata::Dict)
    # Boundaries
    checkopenboundaryphasechange(time, typeinfo, stateinfo, statedata)
    evaluatephasechange(time, typeinfo, stateinfo, statedata)

    # Peridynamics
    evaluatedamage(time, typeinfo, stateinfo, statedata)
end

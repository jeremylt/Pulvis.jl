# ------------------------------------------------------------------------------
# Time step hooks
# ------------------------------------------------------------------------------

"""
```julia
timeloop!(
    particleids,
    typeids,
    neighborindexes,
    dimension,
    bininfo,
    typeinfo,
    stateinfo,
    statedata,
)
```

Evaluate time step; compute density, convert density to pressure, then various forces, acceleration, etc from density and pressure

# Arguments:

  - `particleids::Array{Int}`:     input/output array of local particle ids by type
  - `typeids::Array{Int}`:         input/output array of particle type ids
  - `neighborindexes::Array{Int}`: input/output array of neighbor particle global indexes
  - `dimension::Int`:              spatial dimension of problem
  - `bininfo::BinInfo`:            info about background bin layout
  - `timeinfo::Array{TimeInfo}`:   input/output array with time object
  - `typeinfo::Array{TypeInfo}`:   input/output array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing
"""
function timeloop!(
    particleids::Array{Int},
    typeids::Array{Int},
    neighborindexes::Array{Int},
    bininfo::BinInfo,
    timeinfo::Array{TimeInfo},
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    while continuesimulation()
        # Counter and print header
        preparetimestep!(timeinfo)

        # Pre time step routines
        preevaluatetimestep!(
            bininfo,
            dimension,
            timeinfo[1].time,
            typeinfo,
            stateinfo,
            statedata,
            particleids,
            typeids,
            neighborindexes,
        )

        # Evaluate the physics - the heavy call
        integratetimestep!(
            particleids,
            typeids,
            neighborindexes,
            dimension,
            typeinfo,
            stateinfo,
            statedata,
        )

        # Post time step routines
        postevaluatetimestep!(
            particleids,
            typeids,
            neighborindexes,
            dimension,
            typeinfo,
            stateinfo,
            statedata,
        )

        # Write output and restart
        writeoutput!(typeinfo, stateinfo, statedata)
        # writerestart() # TODO
    end
    return nothing
end

"""
```julia
preparetimestep(timeinfo)
```

Update current time and step number for current time step

# Arguments:

  - `timeinfo::Array{TimeInfo}`: input/output array with time object

# Returns:

  - nothing
"""
function preparetimestep!(timeinfo::Array{TimeInfo})
    timeinfo[1] =
        TimeInfo(timeinfo[1].time + timeinfo[1].dt, timeinfo[1].dt, timeinfo[1].step + 1)
    return nothing
end

"""
```julia
integratetimestep!(
    particleids,
    typeids,
    neighborindexes,
    dimension,
    typeinfo,
    stateinfo,
    statedata,
)
```

Integrate the time step, using the selected interaction model and performing any pre/post routines required by the timestepper

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
function integratetimestep_sph!(
    particleids::Array{Int},
    typeids::Array{Int},
    neighborindexes::Array{Int},
    dimension::Int,
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    # Timestepper pre step routine - TODO

    # Note - updateconnectivitiy method from 4C omitted here because this code is MPI local or run without MPI

    # Interaction model
    evaluatetimestep!(
        particleids,
        typeids,
        neighborindexes,
        dimension,
        typeinfo,
        stateinfo,
        statedata,
    )

    # Timestepper post step routine - TODO
    return nothing
end

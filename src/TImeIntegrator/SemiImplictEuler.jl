# ------------------------------------------------------------------------------
# Time integrator - Semi-implicit Euler
# ------------------------------------------------------------------------------

"""
```julia
preinteractionroutine!(timeinfo, typeinfo, stateinfo, statedata)
```

Time stepper pre particle interaction activities, including advancing the position

# Arguments:

  - `timeinfo::Array{TimeInfo}`:   input/output array with time object
  - `typeinfo::Array{TypeInfo}`:   input/output array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing
"""
function preinteractionroutine_semiimpliciteuler!(timeinfo, typeinfo, stateinfo, statedata)
    hasangularvelocity = findfirst(s -> s.name == "angular velocity", stateinfo) != nothing
    hasangularacceleration =
        findfirst(s -> s.name == "angular acceleration", stateinfo) != nothing

    for i ∈ 1:length(typeinfo)
        # Update velocity and position, then clear acceleration
        if typeinfo[i].name != "open boundary"
            statedata["velocity"][i] .+= timeinfo.dt * statedata["acceleration"][i]
        end
        statedata["position"][i] .+= timeinfo.dt * statedata["velocity"][i]
        statedata["acceleration"[i]] .= 0.0

        # Update angular velocity and acceleration
        if hasangularvelocity && hasangularacceleration
            if length(statedata["angular velocity"][i]) > 0 &&
               length(statedata["angular acceleration"][i]) > 0
                statedata["angular velocity"][i] .+=
                    timeinfo.dt * statedata["angular acceleration"][i]
                statedata["angular acceleration"][i] .= 0.0
            end
        end
    end
    return nothing
end
"""
```julia
postinteractionroutine!(timeinfo, typeinfo, stateinfo, statedata)
```

Time stepper post particle interaction activities

# Arguments:

  - `timeinfo::Array{TimeInfo}`:   input/output array with time object
  - `typeinfo::Array{TypeInfo}`:   input/output array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing
"""
function postinteractionroutine_semiimpliciteuler!(timeinfo, typeinfo, stateinfo, statedata)
    return nothing
end

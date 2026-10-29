# ------------------------------------------------------------------------------
# Drivers
# ------------------------------------------------------------------------------

# source files
include("Drivers/CUDA.jl")
include("Drivers/CPU.jl")

"""
```julia
setdevicecomputation(flag)
```

Set use of device computation

# Arguments:

  - `flag::Bool`: boolean flag to use device computation

# Returns:

  - none

# Example:

```jldoctest
# set device computation
setdevicecomputation(true)

# verify
@assert computedensity! == computedensity_cuda!

# set host computation
setdevicecomputation(false)

# verify
@assert computedensity! == computedensity_cpu!

# output

```
"""
function setdevicecomputation(flag::Bool)
    if flag
        global computedensity! = computedensity_cuda!
        global computepressure! = computepressure_cuda!
        global computeforces! = computeforces_cuda!
    else
        global computedensity! = computedensity_cpu!
        global computepressure! = computepressure_cpu!
        global computeforces! = computeforces_cpu!
    end
    return nothing
end

# default driver - cuda
setdevicecomputation(true)

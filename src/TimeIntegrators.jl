# ------------------------------------------------------------------------------
# Time integrators
# ------------------------------------------------------------------------------

# source files
include("TimeIntegrators/SemiImplicitEuler.jl")

"""
```julia
settimeintegrator(integratortype)
```

Set the time integrator

# Arguments:

  - `integratortype::TimeIntegrator.Type`: time integrator to use

# Returns:

  - nothing

# Example:

```jldoctest
# set time integrator
Pulvis.settimeintegrator(TimeIntegrator.SemiImplictEuler)

# verify
@assert Pulvis.preinteractionroutine! == Pulvis.preinteractionroutine_semiimpliciteuler!;

# output

```
"""
function settimeintegrator(integratortype::TimeIntegrator.Type)
    if integratortype == TimeIntegrator.SemiImplicitEuler
        global preinteractionroutine! = preinteractionroutine_semiimplicteuler!
        global postinteractionroutine! = postinteractionroutine_semiimplicteuler!
    end
    return nothing
end

# Default time integrator - semi-implicit Euler
settimeintegrator(TimeIntegrator.SemiImplicitEuler)

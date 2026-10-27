# ------------------------------------------------------------------------------
# Kernel shape functions
# ------------------------------------------------------------------------------

# source files
include("Kernels/Common.jl")
include("Kernels/CubicSpline.jl")
include("Kernels/QuinticSpline.jl")

"""
```julia
setkernelshape(shape)
```

Set the kernel shape function

# Arguments:

  - `shape::Kernel.Shape`: shape function to use

# Returns:

  - none

# Example:

```jldoctest
# set shape function
setkernelshape(Kernel.CubicSpline)

# verify
@assert w == w_CubicSpline

# set shape function
setkernelshape(Kernel.QuinticSpline)

# verify
@assert w0 == w0_QuinticSpline

# output

```
"""
function setkernelshape(shape::Kernel.Shape)
    if shape == Kernel.CubicSpline
        global smoothinglength = smoothinglength_CubicSpline
        global normalizationconstant = normalizationconstant_CubicSpline
        global w0 = w0_CubicSpline
        global w = w_CubicSpline
        global dwdrij = dwdrij_CubicSpline
        global d2wdrij2 = d2wdrij2_CubicSpline
    elseif shape == Kernel.QuinticSpline
        global smoothinglength = smoothinglength_QuinticSpline
        global normalizationconstant = normalizationconstant_QuinticSpline
        global w0 = w0_QuinticSpline
        global w = w_QuinticSpline
        global dwdrij = dwdrij_QuinticSpline
        global d2wdrij2 = d2wdrij2_QuinticSpline
    end
    return nothing
end

# default kernel - cubic spline
setkernelshape(Kernel.CubicSpline)

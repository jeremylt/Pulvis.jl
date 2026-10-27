# ------------------------------------------------------------------------------
# Kernel shape functions
# ------------------------------------------------------------------------------

# source files
include("Kernel/Common.jl")
include("Kernel/CubicSpline.jl")
include("Kernel/QuinticSpline.jl")

# default kernel - cubic spline
smoothinglength = smoothinglength_CubicSpline;
normalizationconstant = normalizationconstant_CubicSpline;
w0 = w0_CubicSpline;
w = w_CubicSpline;
dwdrij = dwdrij_CubicSpline;
d2wdrij2 = d2wdrij2_CubicSpline;

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
        smoothinglength = smoothinglength_CubicSpline
        normalizationconstant = normalizationconstant_CubicSpline
        w0 = w0_CubicSpline
        w = w_CubicSpline
        dwdrij = dwdrij_CubicSpline
        d2wdrij2 = d2wdrij2_CubicSpline
    elseif shape == Kernel.QuinticSpline
        smoothinglength = smoothinglength_QuinticSpline
        normalizationconstant = normalizationconstant_QuinticSpline
        w0 = w0_QuinticSpline
        w = w_QuinticSpline
        dwdrij = dwdrij_QuinticSpline
        d2wdrij2 = d2wdrij2_QuinticSpline
    end
end

# ------------------------------------------------------------------------------
# kernel shape functions
# ------------------------------------------------------------------------------

# source files
include("kernel/common.jl")
include("kernel/cubic_spline.jl")
include("kernel/quintic_spline.jl")

# default kernel - cubic spline
smoothinglength = smoothinglength_cubicspline;
normalizationconstant = normalizationconstant_cubicspline;
w0 = w0_cubicspline;
w = w_cubicspline;
dwdrij = dwdrij_cubicspline;
d2wdrij2 = d2wdrij2_cubicspline;

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
setkernelshape(Kernel.cubicspline)

# verify
@assert w == w_cubicspline

# set shape function
setkernelshape(Kernel.quinticspline)

# verify
@assert w0 == w0_quinticspline

# output

```
"""
function setkernelshape(shape::Kernel.Shape)
    if shape == Kernel.cubicspline
        smoothinglength = smoothinglength_cubicspline
        normalizationconstant = normalizationconstant_cubicspline
        w0 = w0_cubicspline
        w = w_cubicspline
        dwdrij = dwdrij_cubicspline
        d2wdrij2 = d2wdrij2_cubicspline
    elseif shape == Kernel.quinticspline
        smoothinglength = smoothinglength_quinticspline
        normalizationconstant = normalizationconstant_quinticspline
        w0 = w0_quinticspline
        w = w_quinticspline
        dwdrij = dwdrij_quinticspline
        d2wdrij2 = d2wdrij2_quinticspline
    end
end

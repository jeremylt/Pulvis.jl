# ------------------------------------------------------------------------------
# Kernel functions - common
# ------------------------------------------------------------------------------

"""
```julia
gradwij!(r_ij, support, e_ij, gradW_ij)
```

Compute the gradient of the kernel

# Arguments:

  - `r_ij::Float`:      distance between particles
  - `support::Float`:   support radius
  - `dimension::Float`: spatial dimension
  - `e_ij::Array`:      input array of normalized distance between particles
  - `gradW_ij::Array`:  output array of kernel gradient

# Returns:

  - none

# Example:

```jldoctest
# set cubic kernel
setkernelshape(Kernel.CubicSpline)

# setup
r_ij = 0.2
support = 0.8
h = support / 2.0
e_ij = [0.5, sqrt(3.0) / 2.0, 0.0]
gradW_ij = [0.0, 0.0, 0.0]

# normalization constant
normalization_3d = (1.0 / π) / h^3
q = r_ij / h

# verify
w_unnormalized = (-3.0 * q + 2.25 * q^2) * (1.0 / h)
dimension = 3
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradW_ij ≈ w_unnormalized * normalization_3d * e_ij

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# setup
r_ij = 0.2
support = 0.9
h = support / 3.0
e_ij = [0.5, sqrt(3.0) / 2.0, 0.0]
gradW_ij = [0.0, 0.0, 0.0]

# normalization constant
normalization_3d = 3.0 * (1.0 / π) / (359.0 * h^3)
q = r_ij / h

# verify
w_unnormalized = (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4 - 75.0 * (1.0 - q)^4) * (1.0 / h)
dimension = 3
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradW_ij ≈ w_unnormalized * normalization_3d * e_ij

# output

```
"""
function gradwij!(r_ij, support, dimension::Int, e_ij::Array, gradW_ij::Array)
    gradW_ij .= e_ij * dwdrij(r_ij, support, dimension)
    return nothing
end

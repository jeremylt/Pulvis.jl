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
# setup
r_ij = 1.0
support = 4.1
e_ij = [1.0/3.0, 2.0/3.0, 2.0/3.0]
gradW_ij = [0.0, 0.0, 0.0]

# verify
dimension = 1
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradWij[1] ≈ -0.04907232659332006
dimension = 2
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradWij[1:2] ≈ [-0.016327742187484506, -0.03265548437496901]
dimension = 3
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradWij ≈ [-0.005575326600604468, -0.011150653201208937, -0.011150653201208937]

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# verify
dimension = 1
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradWij[1] ≈ -0.08198535868720497
dimension = 2
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradWij[1:2] ≈ [-0.03355638865445558, -0.06711277730891116]
dimension = 3
gradwij!(r_ij, support, dimension, e_ij, gradW_ij)
@assert gradWij ≈ [-0.014011004628756778, -0.028022009257513556, -0.028022009257513556]

# output

```
"""
function gradwij!(r_ij, support, dimension::Int, e_ij::Array, gradW_ij::Array)
    gradW_ij .= e_ij * dwdrij(r_ij, support, dimension)
    return nothing
end

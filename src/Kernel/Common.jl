# ------------------------------------------------------------------------------
# Kernel functions - common
# ------------------------------------------------------------------------------

"""
```julia
gradwij!(rij, support, eij, gradWij)
```

Compute the gradient of the kernel

# Arguments:

  - `rij::Float`:       distance between particles
  - `support::Float`:   support radius
  - `dimension::Float`: spatial dimension
  - `eij::Array`:       input array of normalized distance between particles
  - `gradWij::Array`:   output array of kernel gradient

# Returns:

  - none

# Example:

```jldoctest
# setup
rij = 1.0
support = 4.1
eij = [1.0/3.0, 2.0/3.0, 2.0/3.0]
gradWij = [0.0, 0.0, 0.0]

# verify
dimension = 1
gradwij!(rij, support, dimension, eij, gradWij)
@assert gradWij[1] ≈ -0.04907232659332006
dimension = 2
gradwij!(rij, support, dimension, eij, gradWij)
@assert gradWij[1:2] ≈ [-0.016327742187484506, -0.03265548437496901]
dimension = 3
gradwij!(rij, support, dimension, eij, gradWij)
@assert gradWij ≈ [-0.005575326600604468, -0.011150653201208937, -0.011150653201208937]

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# verify
dimension = 1
gradwij!(rij, support, dimension, eij, gradWij)
@assert gradWij[1] ≈ -0.08198535868720497
dimension = 2
gradwij!(rij, support, dimension, eij, gradWij)
@assert gradWij[1:2] ≈ [-0.03355638865445558, -0.06711277730891116]
dimension = 3
gradwij!(rij, support, dimension, eij, gradWij)
@assert gradWij ≈ [-0.014011004628756778, -0.028022009257513556, -0.028022009257513556]

# output

```
"""
function gradwij!(rij, support, dimension::Int, eij::Array, gradWij::Array)
    gradWij .= eij * dwdrij(rij, support, dimension)
    return nothing
end

# ------------------------------------------------------------------------------
# common kernel functions
# ------------------------------------------------------------------------------

"""
```julia
gradwij!(rij, support, eij, gradWij)
```

Compute the gradient of the kernel

# Arguments:

  - `rij::Float`:       distance between particles
  - `support::Float`:   support raidus
  - `dimension::Float`: spatial dimension
  - `eij::Array`:       input array of normalized distance between particles
  - `gradWij::Array`:   output array of kernel gradient

# Returns:

  - none

# Example:

```jldoctest
# setup
rij = 3.0
support = 2.1
eij = [1.0/3.0, 2.0/3.0, 2.0/3.0]
gradWij = [0.0, 0.0, 0.0]

# verify
dimension = 1
gradwij!(rij, support, eij, gradWij)
@assert gradwij[1] ≈ 42.0
dimension = 2
gradwij!(rij, support, eij, gradWij)
@assert gradwij[1:2] ≈ [42.0, 42.0]
dimension = 3
gradwij!(rij, support, eij, gradWij)
@assert gradwij ≈ [42.0, 42.0, 42.0]

# set quintic kernel
setkernelshape(Kernel.quinticspline)

# verify
dimension = 1
gradwij!(rij, support, eij, gradWij)
@assert gradwij[1] ≈ 42.0
dimension = 2
gradwij!(rij, support, eij, gradWij)
@assert gradwij[1:2] ≈ [42.0, 42.0]
dimension = 3
gradwij!(rij, support, eij, gradWij)
@assert gradwij ≈ [42.0, 42.0, 42.0]

# output

```
"""
function gradwij!(rij, support, dimension::Int, eij::Array, gradWij::Array)
    gradWij .= eij * dwdrij(rij, support, dimension)
end

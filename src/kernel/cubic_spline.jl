# ------------------------------------------------------------------------------
# cubic spline kernel functions
# ------------------------------------------------------------------------------

πinv = 1.0 / π

"""
```julia
smoothinglength(support)
```

Compute the smoothing length of the kernel

# Arguments:

  - `support::Array`: support raidus

# Returns:

  - smoothing length

# Example:

```jldoctest
# setup
support = 2.1

# set cubic kernel
setkernelshape(Kernel.cubicspline)

# verify
@assert smoothinglength(support) ≈ support / 2.0

# output

```
"""
function smoothinglength_cubicspline(support)
    # support / 2.0
    return 0.5 * support
end

"""
```julia
normalizationconstant(invh, dimension)
```

Compute the normalization constant of the kernel

# Arguments:

  - `invh::Float`:    inverse of the support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - normalization constant

# Example:

```jldoctest
# setup
support = 2.1

# set cubic kernel
setkernelshape(Kernel.cubicspline)

# verify
dimension = 1
@assert normalizationconstant(support, dimension) ≈ (2.0 / 3.0) * (1.0 / support)
dimension = 2
@assert normalizationconstant(support, dimension) ≈ (10.0 / 7.0) / π * (1.0 / support)^2
dimension = 3
@assert normalizationconstant(support, dimension) ≈ (1.0 / π) * (1.0 / support)^3

# output

```
"""
function normalizationconstant_cubicspline(invh, dimension::Int)
    if dimension == 1
        # (2.0 / 3.0) * invh
        return 0.6666666666666666 * invh
    elseif dimension == 2
        # (10.0 / 7.0) / π * invh^2
        return 0.4547284088339866 * invh ^ 2
    elseif dimension == 3
        # (1.0 / π) * invh^3
        return πinv * invh ^ 3
    end
end

"""
```julia
w0(support, dimension)
```

Compute the kernel self interaction

# Arguments:

  - `support::Float`: support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - kernel self interaction

# Example:

```jldoctest
# setup
support = 2.1

# set cubic kernel
setkernelshape(Kernel.cubicspline)

# verify
dimension = 1
@assert w0(support, dimension) ≈ (2.0 / 3.0) * (1.0 / (2.0 / support))
dimension = 2
@assert w0(support, dimension) ≈ (10.0 / 7.0) / π * (1.0 / (2.0 / support))^2
dimension = 3
@assert w0(support, dimension) ≈ (1.0 / π) * (1.0 / (2.0 / support))^3

# output

```
"""
function w0_cubicspline(support, dimension::Int)
    return normalizationconstant_cubicspline(2.0 / support, dimension)
end

"""
```julia
w(rij, support, dimension)
```

Compute the kernel interaction

# Arguments:

  - `rij::Float`:     distance between the points
  - `support::Float`: support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - kernel interaction

# Example:

```jldoctest
# setup
rij     = 3.0::
support = 2.1

# set cubic kernel
setkernelshape(Kernel.cubicspline)

# verify
dimension = 1
@assert w(rij, support, dimension) ≈ 4.2
dimension = 2
@assert w(rij, support, dimension) ≈ 42.0
dimension = 3
@assert w(rij, support, dimension) ≈ 0.0

# output

```
"""
function w_cubicspline(rij, support, dimension::Int)
    invh = 2.0 / support
    q = rij * invh

    if q < 1.0
        return (1.0 - 1.5 * q ^ 2 + 0.75 * q ^ 3) * normalizationconstant_cubicspline(invh, dimension)
    elseif q < 2.0
        return (0.25 * (2.0 - q) ^ 3) * normalizationconstant_cubicspline(invh, dimension)
    else
        return 0.0
    end
end

"""
```julia
dwdrij(rij, support, dimension)
```

Compute the first derivative of the kernel

# Arguments:

  - `rij::Float`:     distance between the points
  - `support::Float`: support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - first derivative of the kernel

# Example:

```jldoctest
# setup
rij     = 3.0
support = 2.1

# set cubic kernel
setkernelshape(Kernel.cubicspline)

# verify
dimension = 1
@assert dwdrij(rij, support, dimension) ≈ 4.2
dimension = 2
@assert dwdrij(rij, support, dimension) ≈ 42.0
dimension = 3
@assert dwdrij(rij, support, dimension) ≈ 0.0

# output

```
"""
function dwdrij_cubicspline(rij, support, dimension::Int)
    invh = 2.0 / support
    q = rij * invh

    if q < 1.0
        return (-3.0 * q + 2.25 * q ^ 2) * invh * normalizationconstant_cubicspline(invh, dimension)
    elseif q < 2.0
        return (-0.75 * (2.0 - q) ^ 2) * invh * normalizationconstant_cubicspline(invh, dimension)
    else
        return 0.0
    end
end

"""
```julia
d2wdrij2(rij, support, dimension)
```

Compute the second derivative of the kernel

# Arguments:

  - `rij::Float`:     distance between the points
  - `support::Float`: support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - second derivative of the kernel

# Example:

```jldoctest
# setup
rij     = 3.0
support = 2.1

# set cubic kernel
setkernelshape(Kernel.cubicspline)

# verify
dimension = 1
@assert d2wdrij2(rij, support, dimension) ≈ 4.2
dimension = 2
@assert d2wdrij2(rij, support, dimension) ≈ 42.0
dimension = 3
@assert d2wdrij2(rij, support, dimension) ≈ 0.0

# output

```
"""
function d2wdrij2_cubicspline(rij, support, dimension)
    invh = 2.0 / support
    q = rij * invh

    if q < 1.0
        return (-3.0 + 4.5 * q) * (invh ^ 2) * normalizationconstant_cubicspline(invh, dimension)
    elseif q < 2.0
        return 1.5 * (2.0 - q) * (invh ^ 2) * normalizationconstant_cubicspline(invh, dimension)
    else
        return 0.0
    end
end

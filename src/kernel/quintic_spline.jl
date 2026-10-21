# ------------------------------------------------------------------------------
# quintic spline kernel functions
# ------------------------------------------------------------------------------

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

# set quintic kernel
setkernelshape(Kernel.quinticspline)

# verify
@assert smoothinglength(support) ≈ support / 2.0

# output

```
"""
function smoothinglength_quinticspline(support)
    # support / 3.0
    return 0.3333333333333333 * support
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

# set quintic kernel
setkernelshape(Kernel.quinticspline)

# verify
dimension = 1
@assert normalizationconstant(support, dimension) ≈ (1.0 / support) / 120.0
dimension = 2
@assert normalizationconstant(support, dimension) ≈ (7.0 / 478.0) / π * (1.0 / support)^2
dimension = 3
@assert normalizationconstant(support, dimension) ≈ (3.0 / 359.0) π * (1.0 / support)^3

# output

```
"""
function normalizationconstant_quinticspline(invh, dimension)
    if dimension == 1
        # invh / 120.0
        return 0.0083333333333333 * invh
    elseif dimension == 2
        # (7.0 / 478.0) / π * invh^2
        return 0.0046614418478797 * invh ^ 2
    elseif dimension == 3
        # (3.0 / 359.0) / π * invh^3
        return 0.0026599711937364 * invh ^ 3
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

# set quintic kernel
setkernelshape(Kernel.quinticspline)

# verify
dimension = 1
@assert w0(support, dimension) ≈ 66.0 * (1.0 / (3.0 / support)) / 120.0
dimension = 2
@assert w0(support, dimension) ≈ 66.0 * (7.0 / 478.0) / π * (1.0 / (3.0 / support))^2
dimension = 3
@assert w0(support, dimension) ≈ 66.0 * (3.0 / 359.0) π * (1.0 / (3.0 / support))^3

# output

```
"""
function w0_quinticspline(support, dimension)
    return 66.0 * normalizationconstant_quinticspline(3.0 / support, dimension)
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
rij = 3.0::support = 2.1

# set quintic kernel
setkernelshape(Kernel.quinticspline)

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
function w_quinticspline(rij, support, dimension)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5 + 15.0 * (1.0 - q)^5) *
               normalizationconstant_quinticspline(invh, dimension)
    elseif q < 2.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5) *
               normalizationconstant_quinticspline(invh, dimension)
    elseif q < 3.0
        return (3.0 - q)^5 * normalizationconstant_quinticspline(invh, dimension)
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
rij = 3.0
support = 2.1

# set quintic kernel
setkernelshape(Kernel.quinticspline)

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
function dwdrij_quinticspline(rij, support, dimension)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4 - 75.0 * (1.0 - q)^4) *
               invh *
               normalizationconstant_quinticspline(invh, dimension)
    elseif q < 2.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4) *
               invh *
               normalizationconstant_quinticspline(invh, dimension)
    elseif q < 3.0
        return -5.0 *
               (3.0 - q)^4 *
               invh *
               normalizationconstant_quinticspline(invh, dimension)
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
rij = 3.0
support = 2.1

# set quintic kernel
setkernelshape(Kernel.quinticspline)

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
function d2wdrij2_quinticspline(rij, support, dimension)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3 + 300.0 * (1.0 - q)^3) *
               (invh ^ 2) *
               normalizationconstant_quinticspline(invh, dimension)
    elseif q < 2.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3) *
               (invh ^ 2) *
               normalizationconstant_quinticspline(invh, dimension)
    elseif q < 3.0
        return 20.0 *
               (3.0 - q)^3 *
               (invh ^ 2) *
               normalizationconstant_quinticspline(invh, dimension)
    else
        return 0.0
    end
end

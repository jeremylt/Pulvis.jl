# ------------------------------------------------------------------------------
# Kernel functions - quintic spline
# ------------------------------------------------------------------------------

"""
```julia
smoothinglength(support)
```

Compute the smoothing length of the kernel

# Arguments:

  - `support::Array`: support radius

# Returns:

  - smoothing length

# Example:

```jldoctest
# setup
support = 2.1

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# verify
@assert smoothinglength(support) ≈ support / 3.0

# output

```
"""
function smoothinglength_QuinticSpline(support)
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
invh = 1.0 / support

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# verify
dimension = 1
@assert normalizationconstant(invh, dimension) ≈ (1.0 / support) / 120.0
dimension = 2
@assert normalizationconstant(invh, dimension) ≈ (7.0 / 478.0) / π * (1.0 / support)^2
dimension = 3
@assert normalizationconstant(invh, dimension) ≈ (3.0 / 359.0) / π * (1.0 / support)^3

# output

```
"""
function normalizationconstant_QuinticSpline(invh, dimension)
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
setkernelshape(Kernel.QuinticSpline)

# verify
dimension = 1
@assert w0(support, dimension) ≈ 66.0 * (3.0 / support) / 120.0
dimension = 2
@assert w0(support, dimension) ≈ 66.0 * (7.0 / 478.0) / π * (3.0 / support)^2
dimension = 3
@assert w0(support, dimension) ≈ 66.0 * (3.0 / 359.0) / π * (3.0 / support)^3

# output

```
"""
function w0_QuinticSpline(support, dimension)
    return 66.0 * normalizationconstant_QuinticSpline(3.0 / support, dimension)
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
rij = 2.0

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# verify
dimension = 1
support = 8.1
@assert w(rij, support, dimension) ≈ 0.12308673182572215
dimension = 1
support = 4.1
@assert w(rij, support, dimension) ≈ 0.050604928977178415
dimension = 2
support = 2.1
@assert w(rij, support, dimension) ≈ 5.660228850077883e-7
dimension = 3
support = 1.1
@assert w(rij, support, dimension) ≈ 0.0

# output

```
"""
function w_QuinticSpline(rij, support, dimension)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5 + 15.0 * (1.0 - q)^5) *
               normalizationconstant_QuinticSpline(invh, dimension)
    elseif q < 2.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5) *
               normalizationconstant_QuinticSpline(invh, dimension)
    elseif q < 3.0
        return (3.0 - q)^5 * normalizationconstant_QuinticSpline(invh, dimension)
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
rij = 2.0

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# verify
dimension = 1
support = 1.1
@assert dwdrij(rij, support, dimension) ≈ 0.0
dimension = 1
support = 8.1
@assert dwdrij(rij, support, dimension) ≈ -0.06306498329433108
dimension = 2
support = 4.1
@assert dwdrij(rij, support, dimension) ≈ -0.04635974146374518
dimension = 3
support = 2.1
@assert dwdrij(rij, support, dimension) ≈ -2.3070798053015142e-5

# output

```
"""
function dwdrij_QuinticSpline(rij, support, dimension)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4 - 75.0 * (1.0 - q)^4) *
               invh *
               normalizationconstant_QuinticSpline(invh, dimension)
    elseif q < 2.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4) *
               invh *
               normalizationconstant_QuinticSpline(invh, dimension)
    elseif q < 3.0
        return -5.0 *
               (3.0 - q)^4 *
               invh *
               normalizationconstant_QuinticSpline(invh, dimension)
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
rij = 2.0

# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# verify
dimension = 1
support = 2.1
@assert d2wdrij2(rij, support, dimension) ≈ 0.0014166432920523402
dimension = 1
support = 1.1
@assert d2wdrij2(rij, support, dimension) ≈ 0.0
dimension = 2
support = 8.1
@assert d2wdrij2(rij, support, dimension) ≈ -0.00032958756972976803
dimension = 3
support = 4.1
@assert d2wdrij2(rij, support, dimension) ≈ 0.03013882982464853

# output

```
"""
function d2wdrij2_QuinticSpline(rij, support, dimension)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3 + 300.0 * (1.0 - q)^3) *
               (invh ^ 2) *
               normalizationconstant_QuinticSpline(invh, dimension)
    elseif q < 2.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3) *
               (invh ^ 2) *
               normalizationconstant_QuinticSpline(invh, dimension)
    elseif q < 3.0
        return 20.0 *
               (3.0 - q)^3 *
               (invh ^ 2) *
               normalizationconstant_QuinticSpline(invh, dimension)
    else
        return 0.0
    end
end

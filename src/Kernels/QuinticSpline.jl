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
# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# setup
support = 0.9
h = support / 3.0

# verify
@assert smoothinglength(support) ≈ h

# output

```
"""
function smoothinglength_QuinticSpline(support)
    # support / 3.0
    return 0.3333333333333333 * support
end

"""
```julia
normalizationconstant(h_inv, dimension)
```

Compute the normalization constant of the kernel

# Arguments:

  - `h_inv::Float`:   inverse of the support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - normalization constant

# Example:

```jldoctest
# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# setup
support = 0.9
h = support / 3.0
h_inv = 1.0 / h

# normalization constants
normalization_1d = 1.0 / (120.0 * h)
normalization_2d = 7.0 * (1.0 / π) / (478.0 * h^2)
normalization_3d = 3.0 * (1.0 / π) / (359.0 * h^3)

# verify
dimension = 1
@assert normalizationconstant(h_inv, dimension) ≈ normalization_1d
dimension = 2
@assert normalizationconstant(h_inv, dimension) ≈ normalization_2d
dimension = 3
@assert normalizationconstant(h_inv, dimension) ≈ normalization_3d

# output

```
"""
function normalizationconstant_QuinticSpline(h_inv, dimension)
    if dimension == 1
        # h_inv / 120.0
        return 0.0083333333333333 * h_inv
    elseif dimension == 2
        # (7.0 / 478.0) / π * h_inv^2
        return 0.0046614418478797 * h_inv ^ 2
    elseif dimension == 3
        # (3.0 / 359.0) / π * h_inv^3
        return 0.0026599711937364 * h_inv ^ 3
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
# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# setup
support = 0.9
h = support / 3.0

# normalization constants
normalization_1d = 1.0 / (120.0 * h)
normalization_2d = 7.0 * (1.0 / π) / (478.0 * h^2)
normalization_3d = 3.0 * (1.0 / π) / (359.0 * h^3)

# verify
w_unnormalized = 66.0
dimension = 1
@assert w0(support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w0(support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w0(support, dimension) ≈ w_unnormalized * normalization_3d

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

  - `r_ij::Float`:    distance between the points
  - `support::Float`: support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - kernel interaction

# Example:

```jldoctest
# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# setup
support = 0.9
h = support / 3.0

# normalization constants
normalization_1d = 1.0 / (120.0 * h)
normalization_2d = 7.0 * (1.0 / π) / (478.0 * h^2)
normalization_3d = 3.0 * (1.0 / π) / (359.0 * h^3)

# verify
r_ij = 0.0
q = r_ij / h
w_unnormalized = 66.0
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.2
q = r_ij / h
w_unnormalized = (3.0 - q)^5 - 6.0 * (2.0 - q)^5 + 15.0 * (1.0 - q)^5
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.5
q = r_ij / h
w_unnormalized = (3.0 - q)^5 - 6.0 * (2.0 - q)^5
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.8
q = r_ij / h
w_unnormalized = (3.0 - q)^5
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.9
q = r_ij / h
w_unnormalized = 0.0
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

# output

```
"""
function w_QuinticSpline(r_ij, support, dimension)
    h_inv = 3.0 / support
    q = r_ij * h_inv

    if q < 1.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5 + 15.0 * (1.0 - q)^5) *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    elseif q < 2.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5) *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    elseif q < 3.0
        return (3.0 - q)^5 * normalizationconstant_QuinticSpline(h_inv, dimension)
    else
        return 0.0
    end
end

"""
```julia
dwdrij(r_ij, support, dimension)
```

Compute the first derivative of the kernel

# Arguments:

  - `r_ij::Float`:    distance between the points
  - `support::Float`: support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - first derivative of the kernel

# Example:

```jldoctest
# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# setup
support = 0.9
h = support / 3.0

# normalization constants
normalization_1d = 1.0 / (120.0 * h)
normalization_2d = 7.0 * (1.0 / π) / (478.0 * h^2)
normalization_3d = 3.0 * (1.0 / π) / (359.0 * h^3)

# verify
r_ij = 0.0
q = r_ij / h
w_unnormalized = 0.0
dimension = 1
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.2
q = r_ij / h
w_unnormalized = (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4 - 75.0 * (1.0 - q)^4) * (1.0 / h)
dimension = 1
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.5
q = r_ij / h
w_unnormalized = (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4) * (1.0 / h)
dimension = 1
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.8
q = r_ij / h
w_unnormalized = -5.0 * (3.0 - q)^4 * (1.0 / h)
dimension = 1
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.9
q = r_ij / h
w_unnormalized = 0.0
dimension = 1
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

# output

```
"""
function dwdrij_QuinticSpline(r_ij, support, dimension)
    h_inv = 3.0 / support
    q = r_ij * h_inv

    if q < 1.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4 - 75.0 * (1.0 - q)^4) *
               h_inv *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    elseif q < 2.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4) *
               h_inv *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    elseif q < 3.0
        return -5.0 *
               (3.0 - q)^4 *
               h_inv *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    else
        return 0.0
    end
end

"""
```julia
d2wdrij2(r_ij, support, dimension)
```

Compute the second derivative of the kernel

# Arguments:

  - `r_ij::Float`:    distance between the points
  - `support::Float`: support radius
  - `dimension::Int`: spatial dimension

# Returns:

  - second derivative of the kernel

# Example:

```jldoctest
# set quintic kernel
setkernelshape(Kernel.QuinticSpline)

# setup
support = 0.9
h = support / 3.0

# normalization constants
normalization_1d = 1.0 / (120.0 * h)
normalization_2d = 7.0 * (1.0 / π) / (478.0 * h^2)
normalization_3d = 3.0 * (1.0 / π) / (359.0 * h^3)

# verify
r_ij = 0.0
q = r_ij / h
w_unnormalized = -4000.0 / 3.0
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.2
q = r_ij / h
w_unnormalized =
    (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3 + 300.0 * (1.0 - q)^3) * (1.0 / h^2)
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.5
q = r_ij / h
w_unnormalized = (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3) * (1.0 / h^2)
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.8
q = r_ij / h
w_unnormalized = 20.0 * (3.0 - q)^3 * (1.0 / h^2)
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.9
q = r_ij / h
w_unnormalized = 0.0
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

# output

```
"""
function d2wdrij2_QuinticSpline(r_ij, support, dimension)
    h_inv = 3.0 / support
    q = r_ij * h_inv

    if q < 1.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3 + 300.0 * (1.0 - q)^3) *
               (h_inv ^ 2) *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    elseif q < 2.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3) *
               (h_inv ^ 2) *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    elseif q < 3.0
        return 20.0 *
               (3.0 - q)^3 *
               (h_inv ^ 2) *
               normalizationconstant_QuinticSpline(h_inv, dimension)
    else
        return 0.0
    end
end

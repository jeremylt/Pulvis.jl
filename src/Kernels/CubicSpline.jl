# ------------------------------------------------------------------------------
# Kernel functions - cubic spline
# ------------------------------------------------------------------------------

π_inv = 1.0 / π

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
# set cubic kernel
setkernelshape(Kernel.CubicSpline)

# setup
support = 0.8
h = support / 2.0

# verify
@assert smoothinglength(support) ≈ h

# output

```
"""
function smoothinglength_CubicSpline(support)
    # support / 2.0
    return 0.5 * support
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
# set cubic kernel
setkernelshape(Kernel.CubicSpline)

# setup
support = 0.4
h = support / 2.0
h_inv = 1.0 / h

# normalization constants
normalization_1d = 2.0 / (3.0 * h)
normalization_2d = 10.0 * (1.0 / π) / (7.0 * h^2)
normalization_3d = (1.0 / π) / h^3

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
function normalizationconstant_CubicSpline(h_inv, dimension::Int)
    if dimension == 1
        # (2.0 / 3.0) * h_inv
        return 0.6666666666666666 * h_inv
    elseif dimension == 2
        # (10.0 / 7.0) / π * h_inv^2
        return 0.4547284088339866 * h_inv ^ 2
    elseif dimension == 3
        # (1.0 / π) * h_inv^3
        return π_inv * h_inv ^ 3
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
# set cubic kernel
setkernelshape(Kernel.CubicSpline)

# setup
support = 0.8
h = support / 2.0

# normalization constants
normalization_1d = 2.0 / (3.0 * h)
normalization_2d = 10.0 * (1.0 / π) / (7.0 * h^2)
normalization_3d = (1.0 / π) / h^3

# verify
w_unnormalized = 1.0
dimension = 1
@assert w0(support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w0(support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w0(support, dimension) ≈ w_unnormalized * normalization_3d

# output

```
"""
function w0_CubicSpline(support, dimension::Int)
    return normalizationconstant_CubicSpline(2.0 / support, dimension)
end

"""
```julia
w(r_ij, support, dimension)
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
# set cubic kernel
setkernelshape(Kernel.CubicSpline)

# setup
support = 0.8
h = support / 2.0

# normalization constants
normalization_1d = 2.0 / (3.0 * h)
normalization_2d = 10.0 * (1.0 / π) / (7.0 * h^2)
normalization_3d = (1.0 / π) / h^3

# verify
r_ij = 0.0
q = r_ij / h
w_unnormalized = 1.0
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.2
q = r_ij / h
w_unnormalized = 1.0 - 1.5 * q^2 + 0.75 * q^3
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.6
q = r_ij / h
w_unnormalized = (2.0 - q)^3 / 4.0
dimension = 1
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert w(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.8
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
function w_CubicSpline(r_ij, support, dimension::Int)
    h_inv = 2.0 / support
    q = r_ij * h_inv

    if q < 1.0
        return (1.0 - 1.5 * q ^ 2 + 0.75 * q ^ 3) *
               normalizationconstant_CubicSpline(h_inv, dimension)
    elseif q < 2.0
        return (0.25 * (2.0 - q) ^ 3) * normalizationconstant_CubicSpline(h_inv, dimension)
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
# set cubic kernel
setkernelshape(Kernel.CubicSpline)

# setup
support = 0.8
h = support / 2.0

# normalization constants
normalization_1d = 2.0 / (3.0 * h)
normalization_2d = 10.0 * (1.0 / π) / (7.0 * h^2)
normalization_3d = (1.0 / π) / h^3

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
w_unnormalized = (-3.0 * q + 2.25 * q^2) * (1.0 / h)
dimension = 1
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.6
q = r_ij / h
w_unnormalized = (-0.75 * (2.0 - q)^2) * (1.0 / h)
dimension = 1
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert dwdrij(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.8
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
function dwdrij_CubicSpline(r_ij, support, dimension::Int)
    h_inv = 2.0 / support
    q = r_ij * h_inv

    if q < 1.0
        return (-3.0 * q + 2.25 * q ^ 2) *
               h_inv *
               normalizationconstant_CubicSpline(h_inv, dimension)
    elseif q < 2.0
        return (-0.75 * (2.0 - q) ^ 2) *
               h_inv *
               normalizationconstant_CubicSpline(h_inv, dimension)
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
# set cubic kernel
setkernelshape(Kernel.CubicSpline)

# setup
support = 0.8
h = support / 2.0

# normalization constants
normalization_1d = 2.0 / (3.0 * h)
normalization_2d = 10.0 * (1.0 / π) / (7.0 * h^2)
normalization_3d = (1.0 / π) / h^3

# verify
r_ij = 0.0
q = r_ij / h
w_unnormalized = -3.0 * (1.0 / h^2)
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.2
q = r_ij / h
w_unnormalized = (-3.0 + 4.5 * q) * (1.0 / h^2)
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.6
q = r_ij / h
w_unnormalized = 1.5 * (2.0 - q) * (1.0 / h^2)
dimension = 1
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_1d
dimension = 2
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_2d
dimension = 3
@assert d2wdrij2(r_ij, support, dimension) ≈ w_unnormalized * normalization_3d

r_ij = 0.8
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
function d2wdrij2_CubicSpline(r_ij, support, dimension::Int)
    h_inv = 2.0 / support
    q = r_ij * h_inv

    if q < 1.0
        return (-3.0 + 4.5 * q) *
               (h_inv ^ 2) *
               normalizationconstant_CubicSpline(h_inv, dimension)
    elseif q < 2.0
        return 1.5 *
               (2.0 - q) *
               (h_inv ^ 2) *
               normalizationconstant_CubicSpline(h_inv, dimension)
    else
        return 0.0
    end
end

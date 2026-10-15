πinv = 1.0 / π

function smoothinglength(support)
    # support / 2.0
    return 0.5 * support
end

function normalizationconstant(invh)
    if dimension == 1
        # (2.0 / 3.0) * invh
        return 0.6666666666666666 * invh
    elseif dimension == 2
        # (10.0 / 7.0) / π * invh^2
        return 0.4547284088339866 * invh ^ 2
    elseif dimension == 3
        # (1.0 / π) & invh^3
        return πinv * invh ^ 3
    end
end

function w0(support)
    return normalizationconstant(2.0 / support)
end

function w(rij, support)
    invh = 2.0 / support
    q = rij * invh

    if q < 1.0
        return (1.0 - 1.5 * q ^ 2 + 0.75 * q ^ 3) * normalizationconstant(invh)
    elseif q < 2.0
        return (0.25 * (2.0 - q) ^ 3) * normalizationconstant(invh)
    else
        return 0.0
    end
end

function dwdrij(rij, support)
    invh = 2.0 / support
    q = rij * invh

    if q < 1.0
        return (-3.0 * q + 2.25 * q ^ 2) * invh * normalizationconstant(invh)
    elseif q < 2.0
        return (-0.75 * (2.0 - q) ^ 2) * invh * normalizationconstant(invh)
    else
        return 0.0
    end
end

function d2wdrij2(rij, support)
    invh = 2.0 / support
    q = rij * invh

    if q < 1.0
        return (-3.0 + 4.5 * q) * (invh ^ 2) * normalizationconstant(invh)
    elseif q < 2.0
        return 1.5 * (2.0 - q) * (invh ^ 2) * normalizationconstant(invh)
    else
        return 0.0
    end
end

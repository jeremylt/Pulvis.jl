function smoothinglength(support)
    # support / 3.0
    return 0.3333333333333333 * support
end

function normalizationconstant(invh)
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

function w0(support)
    return 66.0 * normalizationconstant(3.0 / support)
end

function w(rij, support)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5 + 15.0 * (1.0 - q)^5) * normalizationconstant(invh)
    elseif q < 2.0
        return ((3.0 - q)^5 - 6.0 * (2.0 - q)^5) * normalizationconstant(invh)
    elseif
        return (3.0 - q)^5 * normalizationconstant(invh)
    else
        return 0.0
    end
end

function dwdrij(rij, support)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4 - 75.0 * (1.0 - q)^4) * invh * normalizationconstant(invh)
    elseif q < 2.0
        return (-5.0 * (3.0 - q)^4 + 30.0 * (2.0 - q)^4) * invh * normalizationconstant(invh)
    elseif q < 3.0
        return -5.0 * (3.0 - q)^4 * invh * normalizationconstant(invh)
    else
        return 0.0
    end
end

function d2wdrij2(rij, support)
    invh = 3.0 / support
    q = rij * invh

    if q < 1.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3 + 300.0 *(1.0 - q)^3) * (invh ^ 2) * normalizationconstant(invh)
    elseif q < 2.0
        return (20.0 * (3.0 - q)^3 - 120.0 * (2.0 - q)^3) * (invh ^ 2) * normalizationconstant(invh)
    elseif q < 3.0
        return 20.0 * (3.0 - q)^3 * (invh ^ 2) * normalizationconstant(invh)
    else
        return 0.0
    end
end

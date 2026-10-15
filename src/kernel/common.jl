module julia_sph

# Density

function sumweightedmassself!(particle, radius, mass, density)
    density[particle] = w0(radius[particle]) * mass[particle]
end

function sumweightedmassneighbors!(particle, numneighbors, neighbors, dimension, position, radius, mass, density)
    for neighbor = 1:numneighbors
        neighbor = neighbors[neighborid]
        absdensity = distancebetweenparticles(position[particle * dimension], position[neighbor * dimension])
        density[particle] += w(absdist, radius[particle]) * mass[particle];
    end
end

function computedensity!()
    particle = threadIdx().x
    if particle <= length(density)
        sumweightedmassself!(particle, radius, mass, density)
        sumweightedmassneighbors!(particle, numneighbors, neighbors, dimension, position, radius, mass, density)
    end
end

# Pressure

function computepressure!()
    particle = threadIdx().x
    if particle <= length(pressure)
        densitytopressure!(particle, density, initialdensity, view(pressure, particle:particle+1))
    end
end

# Momentum

function momentumequationneighbors!(particle, numneighbors, neighbors, dimension, typeids, position, radius, mass, density, pressure, velocity, acceleration, modifedvelocity, modifiedacceleration)
    for neighborid = 1:numneighbors
        neighbor = neighbors[neighborid]
        particletype = typeids[particle]
        neighbortype = typeids[neighbor]
        # specific coefficient
        dWdr = 
        speccoeff = specificcoefficient(density[particle], density[neighbor], mass[particle], mass[neighbor], dwdr)
        # pressure gradient
        eij = 
        pressuregradient!(density[particle], density[neighbor], pressure[particle], pressure[neighbor], speccoeff, eij, view(acceleration, particle*dimension:(particle+1)*dimension))
        # shear forces
        absdistance = distancebetweenparticles(position[particle * dimension], position[neighbor * dimension])
        shearforces!(density[particle], density[neighbor], view(velocity, particle*dimension:(particle+1)*dimension), view(velocity, particle*dimension:(particle+1)*dimension), dimension + 2, dynamicviscosity[particletype], dynamicviscosity[neighbortype], bulkviscosity[particletype], bulkviscosity[neighbortype], absdistance, speccoeff, eij, view(acceleration, particle*dimension:(particle+1)*dimension))
        # transport velocity
end

function accelerationaddmomentumcontribution!()
    particle = threadIdx().x
    if particle <= length(pressure)
        momentumequationneighbors!(particle, numneighbors, neighbors, dimension, typeids, position, radius, mass, density, pressure, velocity, acceleration, modifedvelocity, modifiedacceleration)
    end    
end

end # module julia_sph

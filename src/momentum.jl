function momentumequationneighbors!(particle, numneighbors, neighbors, typeids, dynamicviscosity, bulkviscosity, position, radius, mass, density, pressure, velocity, acceleration, modifedvelocity, modifiedacceleration)
    for neighborid = 1:numneighbors
        neighbor = neighbors[neighborid]
        particletype = typeids[particle]
        neighbortype = typeids[neighbor]
        rij = view(position, neighbor*dimension:(neighbor+1)*dimension) - view(position, particle*dimension:(particle+1)*dimension)
        absdistance = distancebetweenparticles(view(position, particle*dimension:(particle+1)*dimension), view(position, neighbor*dimension:(neighbor+1)*dimension))

        # specific coefficient
        dWdrij = dwdrij(absdistance, radius[particle])

        speccoeff = specificcoefficient(density[particle], density[neighbor], mass[particle], mass[neighbor], dwdr)

        # pressure gradient
        eij = (-1.0 / absdistance) * rij

        pressuregradient!(density[particle], density[neighbor], pressure[particle], pressure[neighbor], speccoeff, eij, view(acceleration, particle*dimension:(particle+1)*dimension))

        # shear forces
        shearforces!(density[particle], density[neighbor], view(velocity, particle*dimension:(particle+1)*dimension), view(velocity, neighbor*dimension:(neighbor+1)*dimension), dimension + 2, dynamicviscosity[particletype], dynamicviscosity[neighbortype], bulkviscosity[particletype], bulkviscosity[neighbortype], absdistance, speccoeff, eij, view(acceleration, particle*dimension:(particle+1)*dimension))

        # transport velocity
end

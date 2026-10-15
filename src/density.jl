function sumweightedmassself!(particle, radius, mass, density)
    density[particle] = w0(radius[particle]) * mass[particle]
end

function sumweightedmassneighbors!(particle, numneighbors, neighbors, dimension, position, radius, mass, density)
    for neighbor = 1:numneighbors
        neighbor = neighbors[neighborid]
        absdistance = distancebetweenparticles(view(position, particle*dimension:(particle+1)*dimension), view(position, neighbor*dimension:(neighbor+1)*dimension))

        density[particle] += w(absdistance, radius[particle]) * mass[particle];
    end
end

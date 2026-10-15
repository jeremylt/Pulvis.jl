function gradwij!(rij, support, eij, gradWij)
   gradWij .= eij * dwdrij(rij, support)
end

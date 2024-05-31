function nrm = shape_norm(vec, shp)
    nrm = sqrt(vec.'*shp.shm*vec);
end

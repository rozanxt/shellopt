function grd = shape_gradient(drv, shp)
    grd = zeros(shp.spc.ndof, 1);
    grd(shp.free_dof) = shp.shm(shp.free_dof, shp.free_dof)\drv(shp.free_dof);
end

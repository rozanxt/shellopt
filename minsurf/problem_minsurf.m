function [val, dsc, nrm] = problem_minsurf(dsn, shp)
    shp.nrb = setup_nrb(dsn, shp.nrb);
    shp = assemble_shape(shp);
    [val, drv] = objective_volume(shp);
    
    % Riemannian gradient
    grd = shape_gradient(drv, shp);
    nrm = shape_norm(grd, shp);
    
    % Euclidean gradient
    %grd = drv;
    %grd(shp.drlt_dof) = 0;
    %nrm = norm(grd);
    
    dsc = -grd;
    dsc = reshape(dsc, [], 3)';
    dsc = dsc(:);
end

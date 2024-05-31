function [val, dsc, nrm] = problem_koiter(dsn, geo, shp, sim, pen)
    shp.nrb = setup_nrb(dsn, shp.nrb);
    sim.nrb = nrbrefine(shp.nrb, sim.deg, sim.reg, sim.sub);
    
    sim = assemble_koiter(sim);
    shp = assemble_shape(shp);
    
    [obj, drv] = objective_koiter(sim, shp);
    [ceq, deq] = constraint_volume(geo.acv, shp);
    [val, shd] = penalized_objective(obj, drv, ceq, deq, pen);
    
    % Riemannian gradient
    grd = shape_gradient(shd, shp);
    nrm = shape_norm(grd, shp);
    
    % Euclidean gradient
    %grd = shd;
    %grd(shp.drlt_dof) = 0;
    %nrm = norm(grd);
    
    dsc = -grd;
    dsc = reshape(dsc, [], 3)';
    dsc = dsc(:);
end

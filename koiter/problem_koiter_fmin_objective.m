function [val, grd] = problem_koiter_fmin_objective(dsn, shp, sim)
    shp.nrb = setup_nrb(dsn, shp.nrb);
    sim.nrb = nrbrefine(shp.nrb, sim.deg, sim.reg, sim.sub);
    
    sim = assemble_koiter(sim);
    shp = assemble_shape(shp);
    
    [val, drv] = objective_koiter(sim, shp);
    
    %grd = shape_gradient(drv, shp);
    
    grd = drv;
    grd(shp.drlt_dof) = 0;

    grd = reshape(grd, [], 3)';
    grd = grd(:);
end

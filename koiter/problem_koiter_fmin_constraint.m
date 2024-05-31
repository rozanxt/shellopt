function [cin, ceq, din, deq] = problem_koiter_fmin_constraint(dsn, geo, shp, sim)
    shp.nrb = setup_nrb(dsn, shp.nrb);
    sim.nrb = nrbrefine(shp.nrb, sim.deg, sim.reg, sim.sub);
    
    shp = assemble_shape(shp);
    
    [val, drv] = constraint_volume(geo.acv, shp);
    
    grd = drv;
    grd(shp.drlt_dof) = 0;
    
    grd = reshape(grd, [], 3)';
    grd = grd(:);
    
    cin = [];
    ceq = val;
    din = [];
    deq = grd;
end

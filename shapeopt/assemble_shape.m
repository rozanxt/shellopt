function shp = assemble_shape(shp)
    shp.msh = setup_msh(shp.nrb, shp.ngp);
    shp.spc = setup_spc(shp.nrb, shp.msh);
    [shp.free_dof, shp.drlt_dof] = setup_dof(shp.spc, shp.drlt_bnd, shp.drlt_cps);
    
    %shp.shm = op_u_v_tp(shp.spc, shp.spc, shp.msh)+op_gradu_gradv_tp(shp.spc, shp.spc, shp.msh);
    shp.shm = op_u_v_proj_tp(shp.spc, shp.spc, shp.msh)+op_gradu_gradv_proj_tp(shp.spc, shp.spc, shp.msh);
end

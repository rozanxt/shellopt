function sim = assemble_koiter(sim)
    sim.msh = setup_msh(sim.nrb, sim.ngp);
    sim.spc = setup_spc(sim.nrb, sim.msh);
    [sim.free_dof, sim.drlt_dof] = setup_dof(sim.spc, sim.drlt_bnd, sim.drlt_cps);
    
    sim.ksm = op_KL_shells_tp(sim.spc, sim.spc, sim.msh, sim.E, sim.nu, sim.th);
    sim.rhs = op_f_v_tp(sim.spc, sim.msh, sim.f);
    
    sim.sol = zeros(sim.spc.ndof, 1);
    sim.sol(sim.free_dof) = sim.ksm(sim.free_dof, sim.free_dof)\sim.rhs(sim.free_dof);
end

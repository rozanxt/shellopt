function [csm, lnk] = setup_c1c_vector(sim, mdl)
    csm = zeros(sim.ndf);
    lnk = 1:sim.ndf;
    for fc = 1:numel(mdl.fcs)
        ifc = mdl.fcs{fc};
        pca = ifc(1, 1);
        pcb = ifc(2, 1);
        sda = ifc(1, 2);
        sdb = ifc(2, 2);
        
%         msh = msh_eval_boundary_side(sim.msh{pca}, sda);
%         bda = msh_boundary_side_from_interior(sim.msh{pca}, sda);
%         bdb = msh_boundary_side_from_interior(sim.msh{pcb}, sdb);
%         spa = sp_precompute(sim.spc{pca}.constructor(bda), bda, 'value', true, 'gradient', true);
%         spb = sp_precompute(sim.spc{pcb}.constructor(bdb), bdb, 'value', true, 'gradient', true);
%         
%         if mdl.flp(fc)
%             caa = op_gradu_n_gradv_n(spa, spa, msh, 1);
%             cba = op_gradu_n_gradv_n_flip(spa, spb, msh, 1);
%             cab = cba';
%             cbb = op_gradu_n_gradv_n(spb, spb, msh, 1);
%         else
%             caa = op_gradu_n_gradv_n(spa, spa, msh, 1);
%             cba = op_gradu_n_gradv_n(spa, spb, msh, 1);
%             cab = cba';
%             cbb = op_gradu_n_gradv_n(spb, spb, msh, 1);
%         end
%         
%         csm(sim.dfs{pca}, sim.dfs{pca}) = csm(sim.dfs{pca}, sim.dfs{pca})+caa;
%         csm(sim.dfs{pca}, sim.dfs{pcb}) = csm(sim.dfs{pca}, sim.dfs{pcb})-cab;
%         csm(sim.dfs{pcb}, sim.dfs{pca}) = csm(sim.dfs{pcb}, sim.dfs{pca})-cba;
%         csm(sim.dfs{pcb}, sim.dfs{pcb}) = csm(sim.dfs{pcb}, sim.dfs{pcb})+cbb;
        
        csa = reshape(sim.dfs{pca}, [], 2);
        if sda == 1
            csa = csa(1 : 2, :);
            csa = csa(:);
        elseif sda == 2
            csa = csa((end - 1) : end, :);
            csa = csa(:);
        end
        
        csb = reshape(sim.dfs{pcb}, [], 2);
        if sdb == 1
            csb = csb(1 : 2, :);
            csb = csb(:);
        elseif sdb == 2
            csb = csb((end - 1) : end, :);
            csb = csb(:);
        end
        
        csm(csa, csa) = csm(csa, csa) + blkdiag([1, -1; -1, 1], [1, -1; -1, 1]);
        csm(csb, csb) = csm(csb, csb) + blkdiag([1, -1; -1, 1], [1, -1; -1, 1]);
        csm(csa, csb) = csm(csa, csb) + blkdiag([-1, 1; 1, -1], [-1, 1; 1, -1]);
        csm(csb, csa) = csm(csb, csa) + blkdiag([-1, 1; 1, -1], [-1, 1; 1, -1]);
        
        dfa = sim.spc{pca}.boundary(sda).dofs;
        dfb = sim.spc{pcb}.boundary(sdb).dofs;
%         if mdl.flp(fc)
%             dfb = reshape(flipud(reshape(dfb, [], sim.spc{pcb}.ncomp)), 1, []);
%         end
        itf = [sim.dfs{pca}(dfa); sim.dfs{pcb}(dfb)]';
        lnk(itf) = repmat(min(lnk(itf), [], 2), 1, 2);
    end
end

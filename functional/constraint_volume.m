function [ceq, deq] = constraint_volume(acv, shp)
    mev = msh_precompute(shp.msh);
    sev = sp_precompute(shp.spc, mev, 'value', true, 'gradient', true, 'divergence', true);
    
    wjs = mev.quad_weights.*mev.jacdet;
    
    ceq = acv-sum(sum(wjs));
    
    if nargout >= 2
        deq = zeros(shp.spc.ndof, 1);
        for el = 1:mev.nel
            for qn = 1:mev.nqn
                wj = wjs(qn, el);
                for sh = 1:sev.nsh(el)
                    dof = sev.connectivity(sh, el);
                    sfd = sev.shape_function_divs(qn, sh, el);
                    deq(dof) = deq(dof)-wj*sfd;
                end
            end
        end
    end
end

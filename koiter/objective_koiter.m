function [obj, drv] = objective_koiter(sim, shp)
    obj = 0.5*sim.sol.'*sim.ksm*sim.sol;
    
    if nargout >= 2
        mev = msh_precompute(shp.msh);
        wjs = mev.quad_weights.*mev.jacdet;
        
        kev = sp_precompute_param(sim.spc, shp.msh, 'value', true, 'gradient', true, 'hessian', true);
        [sms, ems] = op_KL_membrane_stress(kev, kev, mev, sim.E(0, 0, 0), sim.nu(0, 0, 0), sim.th);
        [sbs, ebs] = op_KL_bending_stress(kev, kev, mev, sim.E(0, 0, 0), sim.nu(0, 0, 0), sim.th);
        
        sev = sp_precompute_param(shp.spc, shp.msh, 'value', true, 'gradient', true, 'hessian', true);
        
        drv = zeros(sev.ndof, 1);
        for el = 1:mev.nel
            for qn = 1:mev.nqn
                wj = wjs(qn, el);
                
                dof = kev.connectivity(:, el);
                uj = sum(repmat(permute(sim.sol(dof), [3, 2, 1]), 3, 2).*squeeze(kev.shape_function_gradients(:, :, qn, :, el)), 3);
                uh = sum(repmat(permute(sim.sol(dof), [4, 2, 3, 1]), 3, 2, 2).*squeeze(kev.shape_function_hessians(:, :, :, qn, :, el)), 4);
                sm = sum(repmat(permute(sim.sol(dof), [3, 2, 1]), 2, 2).*squeeze(sms(:, :, qn, :, el)), 3);
                em = sum(repmat(permute(sim.sol(dof), [3, 2, 1]), 2, 2).*squeeze(ems(:, :, qn, :, el)), 3);
                sb = sum(repmat(permute(sim.sol(dof), [3, 2, 1]), 2, 2).*squeeze(sbs(:, :, qn, :, el)), 3);
                eb = sum(repmat(permute(sim.sol(dof), [3, 2, 1]), 2, 2).*squeeze(ebs(:, :, qn, :, el)), 3);
                
                jac = mev.geo_map_jac(:, :, qn, el);
                sfm = jac.'*jac;
                nrm = mev.normal(:, qn, el);
                axa = cross(jac(:, 1), jac(:, 2));
                len = norm(axa);
                
                uj_arr = repmat(uj, 1, 1, 2);
                jac_arr = repmat(jac, 1, 1, 2);
                nrm_arr = repmat(nrm, 1, 2, 2);
                
                axa_dot_u = cross(uj_arr, permute(jac_arr, [1, 3, 2]))+cross(jac_arr, permute(uj_arr, [1, 3, 2]));
                n_dot_u = (axa_dot_u-dot(nrm_arr, axa_dot_u).*nrm_arr)/len;
                
                for sh = 1:sev.nsh(el)
                    dof = sev.connectivity(sh, el);
                    sfj = sev.shape_function_gradients(:, :, qn, sh, el);
                    sfh = sev.shape_function_hessians(:, :, :, qn, sh, el);
                    
                    sfj_arr = repmat(sfj, 1, 1, 2);
                    
                    em_dot = squeeze(0.5*(dot(uj_arr, permute(sfj_arr, [1, 3, 2]))+dot(sfj_arr, permute(uj_arr, [1, 3, 2]))));
                    axa_dot_sf = cross(sfj_arr, permute(jac_arr, [1, 3, 2]))+cross(jac_arr, permute(sfj_arr, [1, 3, 2]));
                    n_dot_sf = (axa_dot_sf-dot(nrm_arr, axa_dot_sf).*nrm_arr)/len;
                    eb_dot = squeeze(-dot(uh, n_dot_sf)-dot(sfh, n_dot_u));
                    
                    drv(dof) = drv(dof)-wj*(sum(sum(sm.*em_dot+sb.*eb_dot))+0.25*sum(sum(sm.*em+sb.*eb))*trace(sfm\(sfj.'*jac+jac.'*sfj)));
                end
            end
        end
    end
end

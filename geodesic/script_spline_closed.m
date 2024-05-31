A = 1.0;
dt = 0.1;
m = 100;

sim.deg = 3;
sim.reg = 2;
sim.sub = 8;
sim.ngp = 4;

vis.sub = 200;

nrb_pts = zeros(4, 3);
nrb_pts(:, 1) = [1; 0; 0; 1];
nrb_pts(:, 2) = [1; 1; 0; 1] / sqrt(2);
nrb_pts(:, 3) = [0; 1; 0; 1];
nrb_knt = [0, 0, 0, 1, 1, 1];
nrb = nrbmak(nrb_pts, nrb_knt);

mdl.nrb = {nrb,...
           nrbtform(nrb, vecrotz(pi / 2)),...
           nrbtform(nrb, vecrotz(pi)),...
           nrbtform(nrb, vecrotz(3 * pi / 2))};

mdl.fcs = {[1, 2; 2, 1], [2, 2; 3, 1], [3, 2; 4, 1], [4, 2; 1, 1]};

sim.nrb = cell(1, numel(mdl.nrb));
sim.msh = cell(1, numel(mdl.nrb));
spc_scl.spc = cell(1, numel(mdl.nrb));
spc_scl.dfs = cell(1, numel(mdl.nrb));
spc_scl.ndf = 0;
spc_vec.spc = cell(1, numel(mdl.nrb));
spc_vec.dfs = cell(1, numel(mdl.nrb));
spc_vec.ndf = 0;
for pc = 1 : numel(mdl.nrb)
    nrb = nrbrefine(mdl.nrb{pc}, sim.deg, sim.reg, sim.sub);
    msh = setup_msh(nrb, sim.ngp);
    sps = sp_nurbs(nrb, msh);
    spv = sp_vector(repmat({sps}, 1, msh.rdim), msh);
    
    sim.nrb{pc} = nrb;
    sim.msh{pc} = msh;
    
    spc_scl.spc{pc} = sps;
    spc_scl.dfs{pc} = (spc_scl.ndf + 1) : (spc_scl.ndf + sps.ndof);
    spc_scl.ndf = spc_scl.ndf + sps.ndof;
    
    spc_vec.spc{pc} = spv;
    spc_vec.dfs{pc} = (spc_vec.ndf + 1) : (spc_vec.ndf + spv.ndof);
    spc_vec.ndf = spc_vec.ndf + spv.ndof;
end
clear pc nrb msh sps spv;

[spc_scl.csm, spc_scl.lnk] = setup_c1c_scalar(spc_scl, mdl);

spc_scl.dof = unique(spc_scl.lnk);

idm = eye(spc_scl.ndf);
spc_scl.ns0 = idm(spc_scl.lnk, :);
spc_scl.ns0 = spc_scl.ns0(:, spc_scl.dof);
spc_scl.ns1 = null(spc_scl.ns0.' * spc_scl.csm * spc_scl.ns0, 'r');
spc_scl.nsm = spc_scl.ns0 * spc_scl.ns1;

fcn = @(x, y) sin(5 * atan2(y, x));%ones(size(x));
mat = zeros(spc_scl.ndf, spc_scl.ndf);
rhs = zeros(spc_scl.ndf, 1);
for pc = 1 : numel(mdl.nrb)
    msh = sim.msh{pc};
    spc = spc_scl.spc{pc};
    dof = spc_scl.dfs{pc};
    
    mat(dof, dof) = op_u_v_tp(spc, spc, msh);
    rhs(dof) = op_f_v_tp(spc, msh, fcn);
end
clear pc msh spc dof;

mat = spc_scl.nsm.' * mat * spc_scl.nsm;
rhs = spc_scl.nsm.' * rhs;
sol = mat \ rhs;
anf = spc_scl.nsm * sol;

cfg = zeros(spc_vec.ndf, 1);
for pc = 1 : numel(mdl.nrb)
    nrb = sim.nrb{pc};
    dof = spc_vec.dfs{pc};
    
    cfs = permute(nrb.coefs(1 : 2, :, :) ./ nrb.coefs(4, :, :), [2, 3, 1]);
    cfg(dof) = cfs(:);
end
clear pc nrb dof cfs;

figure;
for i = 1 : m
    % Set up multipatch isogeometric space
    sim.nrb = cell(1, numel(mdl.nrb));
    sim.msh = cell(1, numel(mdl.nrb));
    spc_scl.spc = cell(1, numel(mdl.nrb));
    spc_scl.dfs = cell(1, numel(mdl.nrb));
    spc_scl.ndf = 0;
    spc_vec.spc = cell(1, numel(mdl.nrb));
    spc_vec.dfs = cell(1, numel(mdl.nrb));
    spc_vec.ndf = 0;
    for pc = 1 : numel(mdl.nrb)
        nrb = nrbrefine(mdl.nrb{pc}, sim.deg, sim.reg, sim.sub);
        msh = setup_msh(nrb, sim.ngp);
        sps = sp_nurbs(nrb, msh);
        spv = sp_vector(repmat({sps}, 1, msh.rdim), msh);
        
        sim.nrb{pc} = nrb;
        sim.msh{pc} = msh;
        
        spc_scl.spc{pc} = sps;
        spc_scl.dfs{pc} = (spc_scl.ndf + 1) : (spc_scl.ndf + sps.ndof);
        spc_scl.ndf = spc_scl.ndf + sps.ndof;
        
        spc_vec.spc{pc} = spv;
        spc_vec.dfs{pc} = (spc_vec.ndf + 1) : (spc_vec.ndf + spv.ndof);
        spc_vec.ndf = spc_vec.ndf + spv.ndof;
    end
    clear pc nrb msh sps spv;
    
    [spc_scl.csm, spc_scl.lnk] = setup_c1c_scalar(spc_scl, mdl);
    
    spc_scl.dof = unique(spc_scl.lnk);
    
    idm = eye(spc_scl.ndf);
    spc_scl.ns0 = idm(spc_scl.lnk, :);
    spc_scl.ns0 = spc_scl.ns0(:, spc_scl.dof);
    spc_scl.ns1 = null(spc_scl.ns0.' * spc_scl.csm * spc_scl.ns0, 'r');
    spc_scl.nsm = spc_scl.ns0 * spc_scl.ns1;
    
    [spc_vec.csm, spc_vec.lnk] = setup_c1c_vector(spc_vec, mdl);
    
    spc_vec.dof = unique(spc_vec.lnk);
    
    idm = eye(spc_vec.ndf);
    spc_vec.ns0 = idm(spc_vec.lnk, :);
    spc_vec.ns0 = spc_vec.ns0(:, spc_vec.dof);
    spc_vec.ns1 = null(spc_vec.ns0.' * spc_vec.csm * spc_vec.ns0, 'r');
    spc_vec.nsm = spc_vec.ns0 * spc_vec.ns1;
    
    % Compute velocity vector field
    mat = zeros(spc_vec.ndf, spc_vec.ndf);
    rhs = zeros(spc_vec.ndf, 1);
    for pc = 1 : numel(mdl.nrb)
        msh = sim.msh{pc};
        sps = spc_scl.spc{pc};
        spv = spc_vec.spc{pc};
        dof = spc_vec.dfs{pc};
        
        mev = msh_precompute(msh);
        ses = sp_precompute_param(sps, msh);
        sev = sp_precompute_param(spv, msh);
        
        jac = squash(mev.geo_map_jac);
        nrm = [jac(2, :, :); -jac(1, :, :)] ./ permute(mev.jacdet, [3, 1, 2]);
        
        mat0 = op_u_v_tp(spv, spv, msh);
        mat1 = op_gradu_gradv_tp(spv, spv, msh);
        
        anv = squash(sum(permute(anf(ses.connectivity + (pc - 1) * ses.ndof), [3, 1, 2]) .* ses.shape_functions, 2));
        anb = op_f_v(sev, mev, permute(anv, [3, 1, 2]) .* nrm);
        
        mat(dof, dof) = mat0 + A * mat1;
        rhs(dof) = anb;
    end
    clear pc msh sps spv dof mev ses sev jac nrm mat0 mat1 anv anb;
    
    mat = spc_vec.nsm.' * mat * spc_vec.nsm;
    rhs = spc_vec.nsm.' * rhs;
    sol = mat \ rhs;
    vel = spc_vec.nsm * sol;
    
    % Compute change in momentum
    mat2 = zeros(spc_scl.ndf, spc_scl.ndf);
    rhs2 = zeros(spc_scl.ndf, 1);
    for pc = 1 : numel(mdl.nrb)
        msh = sim.msh{pc};
        sps = spc_scl.spc{pc};
        spv = spc_vec.spc{pc};
        dof = spc_scl.dfs{pc};
        
        mev = msh_precompute(msh);
        ses = sp_precompute_param(sps, msh);
        sev = sp_precompute_param(spv, msh, 'value', true, 'gradient', true);
        
        jac = squash(mev.geo_map_jac);
        hes = squash(mev.geo_map_der2);
        nrm = [jac(2, :, :); -jac(1, :, :)] ./ permute(mev.jacdet, [3, 1, 2]);
        jc2 = squash(dot(jac, jac));
        kpc = -squash(dot(nrm, hes)) ./ jc2;
        
        vef = permute(vel(sev.connectivity + (pc - 1) * sev.ndof), [3, 4, 1, 2]);
        vct = squash(sum(vef .* sev.shape_functions, 3));
        ct2 = squash(dot(vct, vct));
        
        veg = permute(vel(sev.connectivity + (pc - 1) * sev.ndof), [3, 4, 5, 1, 2]);
        dct = squash(sum(veg .* sev.shape_function_gradients, 4));
        dsct2 = squash(dot(dct, dct)) ./ jc2;
        
        bnv = (0.5 * kpc) .* (A * dsct2 - ct2);
        
        mat2(dof, dof) = op_u_v_tp(sps, sps, msh);
        rhs2(dof) = op_f_v(ses, mev, bnv);
    end
    clear pc msh sps spv dof mev ses sev jac hes nrm jc2 kpc vef vct ct2 veg dct dsct2 trj anv bnv;
    
    mat2 = spc_scl.nsm.' * mat2 * spc_scl.nsm;
    rhs2 = spc_scl.nsm.' * rhs2;
    sol2 = mat2 \ rhs2;
    bnf = spc_scl.nsm * sol2;
    
    cfg = cfg + vel * dt;
    anf = anf + bnf * dt;
    
    for pc = 1 : numel(mdl.nrb)
        sim.nrb{pc}.coefs(1 : 2, :, :) = reshape(cfg(spc_vec.dfs{pc}), size(sim.nrb{pc}.coefs, 2), sim.msh{pc}.rdim).' .* sim.nrb{pc}.coefs(4, :, :);
    end
    
%     me0 = msh_precompute(setup_msh(nrb0, sim.ngp));
%     tst = reshape(permute(reshape(bnv, sim.ngp(1), sim.ngp(2), sim.sub(1), sim.sub(2)), [1, 3, 2, 4]), sim.ngp(1) * sim.sub(1), sim.ngp(2) * sim.sub(2));
%     mst = reshape(permute(reshape(me0.geo_map, 3, sim.ngp(1), sim.ngp(2), sim.sub(1), sim.sub(2)), [1, 2, 4, 3, 5]), 3, sim.ngp(1) * sim.sub(1), sim.ngp(2) * sim.sub(2));
%     [X, Y] = deal(squash(mst(1, :, :)), squash(mst(2, :, :)));
%     surf(X, Y, tst);
    
%     [val, grd] = sp_eval(bnf, spc_scalar, geo_load(nrb0), vis.sub);
%     [X, Y] = deal(squash(grd(1, :, :)), squash(grd(2, :, :)));
%     surf(X, Y, val);
    
    for pc = 1 : numel(mdl.nrb)
        nrbplot(sim.nrb{pc}, vis.sub);
        hold on;
    end
    hold off;
    axis equal;
    xlim([-3, 3]);
    ylim([-3, 3]);
    pause;
end

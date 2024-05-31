A = 1.0;
dt = 0.05;
m = 100;

sim.deg = 3;
sim.reg = 2;
sim.sub = 8;
sim.ngp = 4;

vis.sub = 200;

nrb_pts = zeros(4, 2);
nrb_pts(:, 1) = [0, 1, 0, 1]';
nrb_pts(:, 2) = [1, 1, 0, 1]';
nrb_knt = [0, 0, 1, 1];
nrb0 = nrbmak(nrb_pts, nrb_knt);
nrb0 = nrbrefine(nrb0, sim.deg, sim.reg, sim.sub);
nrb = nrb0;

msh = setup_msh(nrb, sim.ngp);

spc_scalar = sp_nurbs(nrb, msh);
spc_vector = sp_vector(repmat({spc_scalar}, 1, msh.rdim), msh);

fcn = @(x, y) sin(pi * x);
mat = op_u_v_tp(spc_scalar, spc_scalar, msh);
rhs = op_f_v_tp(spc_scalar, msh, fcn);
anf = mat \ rhs;

cfg = permute(nrb.coefs(1 : 2, :, :), [2, 3, 1]);
cfg = cfg(:);

figure;
for i = 1 : m
    msh = setup_msh(nrb, sim.ngp);
    mev = msh_precompute(msh);
    
    jac = squeeze(mev.geo_map_jac);
    hes = squeeze(mev.geo_map_der2);
    jcd = mev.jacdet;
    jc2 = squeeze(dot(jac, jac));
    nrm = [-jac(2, :, :); jac(1, :, :)] ./ permute(jcd, [3, 1, 2]);
    kpc = squeeze(dot(nrm, hes)) ./ jc2;
    
    spc_scalar = sp_nurbs(nrb, msh);
    spc_vector = sp_vector(repmat({spc_scalar}, 1, msh.rdim), msh);
    
    sev_scalar = sp_precompute_param(spc_scalar, msh);
    sev_vector = sp_precompute_param(spc_vector, msh, 'value', true, 'gradient', true);
    
    mat0 = op_u_v_tp(spc_vector, spc_vector, msh);
    mat1 = op_gradu_gradv_tp(spc_vector, spc_vector, msh);
    
    anv = squeeze(sum(permute(anf(sev_scalar.connectivity), [3, 1, 2]) .* sev_scalar.shape_functions, 2));
    anb = op_f_v(sev_vector, mev, permute(anv, [3, 1, 2]) .* nrm);
    
    dof = reshape(1 : spc_vector.ndof, [], 2);
    dof = dof(2 : (end - 1), :);
    dof = dof(:);
    sol = (mat0(dof, dof) + A * mat1(dof, dof)) \ anb(dof);
    vel = zeros(spc_vector.ndof, 1);
    vel(dof) = sol;
    
    disp(vel.' * (mat0 + A * mat1) * vel);
    
    vef = permute(vel(sev_vector.connectivity), [3, 4, 1, 2]);
    vct = squeeze(sum(vef .* sev_vector.shape_functions, 3));
    ct2 = squeeze(dot(vct, vct));
    
    veg = permute(vel(sev_vector.connectivity), [3, 4, 5, 1, 2]);
    dct = squeeze(sum(veg .* sev_vector.shape_function_gradients, 4));
    dsct2 = squeeze(dot(dct, dct)) ./ jc2;
    
    trj = squeeze(dot(jac, dct)) ./ jc2;
    
    bnv = (0.5 * kpc) .* (A * dsct2 - ct2) - anv .* trj;
    
    mat2 = op_u_v_tp(spc_scalar, spc_scalar, msh);
    rhs2 = op_f_v(sev_scalar, mev, bnv);
    bnf = mat2 \ rhs2;
    
    cfg = cfg + vel * dt;
    anf = anf + bnf * dt;
    
    cfs = permute(reshape(cfg, size(nrb.coefs, 2), size(nrb.coefs, 3), mev.rdim), [3, 1, 2]);
    nrb.coefs(1 : 2, :, :) = cfs;
    
%     me0 = msh_precompute(setup_msh(nrb0, sim.ngp));
%     tst = reshape(permute(reshape(bnv, sim.ngp(1), sim.ngp(2), sim.sub(1), sim.sub(2)), [1, 3, 2, 4]), sim.ngp(1) * sim.sub(1), sim.ngp(2) * sim.sub(2));
%     mst = reshape(permute(reshape(me0.geo_map, 3, sim.ngp(1), sim.ngp(2), sim.sub(1), sim.sub(2)), [1, 2, 4, 3, 5]), 3, sim.ngp(1) * sim.sub(1), sim.ngp(2) * sim.sub(2));
%     [X, Y] = deal(squeeze(mst(1, :, :)), squeeze(mst(2, :, :)));
%     surf(X, Y, tst);
    
%     [val, grd] = sp_eval(bnf, spc_scalar, geo_load(nrb0), vis.sub);
%     [X, Y] = deal(squeeze(grd(1, :, :)), squeeze(grd(2, :, :)));
%     surf(X, Y, val);
    
    nrbplot(nrb, vis.sub);
    axis equal;
    pause;
end

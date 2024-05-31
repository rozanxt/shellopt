A = 1.0;
dt = 0.05;
m = 100;

sim.deg = [3, 3];
sim.reg = [2, 2];
sim.sub = [8, 8];
sim.ngp = [4, 4];

vis.sub = [20, 20];

nrb0 = nrbtform(unit_square, vecscale([pi, pi, 1]));
nrb0 = nrbrefine(nrb0, sim.deg, sim.reg, sim.sub);
nrb = nrb0;

msh = setup_msh(nrb, sim.ngp);

spc_scalar = sp_nurbs(nrb, msh);
spc_vector = sp_vector(repmat({spc_scalar}, 1, msh.rdim), msh);

fcn = @(x, y, z) sin(x) .* sin(y);
mat = op_u_v_tp(spc_scalar, spc_scalar, msh);
rhs = op_f_v_tp(spc_scalar, msh, fcn);
anf = mat \ rhs;

cfg = permute(nrb.coefs(1 : 3, :, :), [2, 3, 1]);
cfg = cfg(:);

figure;
for i = 1 : m
    msh = setup_msh(nrb, sim.ngp);
    mev = msh_precompute(msh);
    
    jac = mev.geo_map_jac;
    hes = mev.geo_map_der2;
    nrm = mev.normal;
    
    f11 = dot(jac(:, 1, :, :), jac(:, 1, :, :));
    f22 = dot(jac(:, 2, :, :), jac(:, 2, :, :));
    f12 = dot(jac(:, 1, :, :), jac(:, 2, :, :));
    f21 = dot(jac(:, 2, :, :), jac(:, 1, :, :));
    fff = [f11, f12; f21, f22];
    
    fff_inv = [fff(2, 2, :, :), -fff(1, 2, :, :); -fff(2, 1, :, :), fff(1, 1, :, :)];
    fff_inv = fff_inv ./ (fff(1, 1, :, :) .* fff(2, 2, :, :) - fff(1, 2, :, :) .* fff(2, 1, :, :));
    
    sff = -squeeze(dot(hes, repmat(permute(nrm, [1, 4, 5, 2, 3]), 1, 2, 2, 1, 1)));
    
    ais = multiprod(fff_inv, sff);
    trs = squeeze(ais(1, 1, :, :) + ais(2, 2, :, :));
    
    spc_scalar = sp_nurbs(nrb, msh);
    spc_vector = sp_vector(repmat({spc_scalar}, 1, msh.rdim), msh);
    
    sev_scalar = sp_precompute_param(spc_scalar, msh);
    sev_vector = sp_precompute_param(spc_vector, msh, 'value', true, 'gradient', true);
    
    mat0 = op_u_v_tp(spc_vector, spc_vector, msh);
    mat1 = op_gradu_gradv_tp(spc_vector, spc_vector, msh);
    
    anv = squeeze(sum(permute(anf(sev_scalar.connectivity), [3, 1, 2]) .* sev_scalar.shape_functions, 2));
    anb = op_f_v(sev_vector, mev, permute(anv, [3, 1, 2]) .* nrm);
    
    dof = setup_dof(spc_vector, [1, 2, 3, 4], {[1, 2, 3], [1, 2, 3], [1, 2, 3], [1, 2, 3]});
    sol = (mat0(dof, dof) + A * mat1(dof, dof)) \ anb(dof);
    vel = zeros(spc_vector.ndof, 1);
    vel(dof) = sol;
    
    vef = permute(vel(sev_vector.connectivity), [3, 4, 1, 2]);
    vft = squeeze(sum(vef .* sev_vector.shape_functions, 3));
    ft2 = squeeze(dot(vft, vft));
    
    veg = permute(vel(sev_vector.connectivity), [3, 4, 5, 1, 2]);
    gft = squeeze(sum(veg .* sev_vector.shape_function_gradients, 4));
    
    g11 = dot(gft(:, 1, :, :), gft(:, 1, :, :));
    g22 = dot(gft(:, 2, :, :), gft(:, 2, :, :));
    g12 = dot(gft(:, 1, :, :), gft(:, 2, :, :));
    g21 = dot(gft(:, 2, :, :), gft(:, 1, :, :));
    gff = [g11, g12; g21, g22];
    
    aig = multiprod(fff_inv, gff);
    trg = squeeze(aig(1, 1, :, :) + aig(2, 2, :, :));
    
    gsg = multiprod(ais, aig);
    tsg = squeeze(gsg(1, 1, :, :) + gsg(2, 2, :, :));
    
    j11 = dot(jac(:, 1, :, :), gft(:, 1, :, :));
    j22 = dot(jac(:, 2, :, :), gft(:, 2, :, :));
    j12 = dot(jac(:, 1, :, :), gft(:, 2, :, :));
    j21 = dot(jac(:, 2, :, :), gft(:, 1, :, :));
    jff = [j11, j12; j21, j22];
    
    aij = multiprod(fff_inv, jff);
    trj = squeeze(aij(1, 1, :, :) + aij(2, 2, :, :));
    
    bnv = A * tsg - (0.5 * trs) .* (ft2 + A * trg) - anv .* trj;
    
    mat2 = op_u_v_tp(spc_scalar, spc_scalar, msh);
    rhs2 = op_f_v(sev_scalar, mev, bnv);
    bnf = mat2 \ rhs2;
    
    cfg = cfg + vel * dt;
    anf = anf + bnf * dt;
    
    cfs = permute(reshape(cfg, size(nrb.coefs, 2), size(nrb.coefs, 3), mev.rdim), [3, 1, 2]);
    nrb.coefs(1 : 3, :, :) = cfs;
    
%     me0 = msh_precompute(setup_msh(nrb0, sim.ngp));
%     tst = reshape(permute(reshape(bnv, sim.ngp(1), sim.ngp(2), sim.sub(1), sim.sub(2)), [1, 3, 2, 4]), sim.ngp(1) * sim.sub(1), sim.ngp(2) * sim.sub(2));
%     mst = reshape(permute(reshape(me0.geo_map, 3, sim.ngp(1), sim.ngp(2), sim.sub(1), sim.sub(2)), [1, 2, 4, 3, 5]), 3, sim.ngp(1) * sim.sub(1), sim.ngp(2) * sim.sub(2));
%     [X, Y] = deal(squeeze(mst(1, :, :)), squeeze(mst(2, :, :)));
%     surf(X, Y, tst);
    
%     [val, grd] = sp_eval(bnf, spc_scalar, geo_load(nrb0), vis.sub);
%     [X, Y] = deal(squeeze(grd(1, :, :)), squeeze(grd(2, :, :)));
%     surf(X, Y, val);
    
    nrbplot(nrb, vis.sub);
    hold on;
    quiver3(linspace(0, pi, vis.sub(1) + 1).', pi * ones(vis.sub(1) + 1, 1), ones(vis.sub(1) + 1, 1), zeros(vis.sub(1) + 1, 1), zeros(vis.sub(1) + 1, 1), sin(linspace(0, pi, vis.sub(1) + 1).'), 2, 'Color', [0.72549, 0.15686, 0.09804]);
    quiver3(pi * ones(vis.sub(2) + 1, 1), linspace(0, pi, vis.sub(2) + 1).', ones(vis.sub(2) + 1, 1), zeros(vis.sub(2) + 1, 1), zeros(vis.sub(2) + 1, 1), sin(linspace(0, pi, vis.sub(2) + 1).'), 2, 'Color', [0.72549, 0.15686, 0.09804]);
    hold off;
    axis equal;
    xlim([0, pi]);
    ylim([0, pi]);
    zlim([1, pi]);
    %if ismember(i, 4 : 4 : m)
    %    exportgraphics(gcf, ['surface-', int2str(i / 4), '.png'], 'Resolution', 300);
    %end
    pause;
end

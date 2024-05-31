%% Shape optimization of a thin elastic shell using MATLAB's fmincon

% Problem parameters
E = 7e5;
nu = 0.33;
th = 0.1;
f = 5e-2;

fx = @(x, y, z) zeros(size(x));
fy = @(x, y, z) zeros(size(x));
fz = @(x, y, z) -f*ones(size(x));

sim.E = @(x, y, z) E*ones(size(x));
sim.nu = @(x, y, z) nu*ones(size(x));
sim.th = th;
sim.f = @(x, y, z, ind) cat(1,...
    reshape(fx(x, y, z), [1, size(x)]),...
    reshape(fy(x, y, z), [1, size(x)]),...
    reshape(fz(x, y, z), [1, size(x)]));

sim.drlt_bnd = [1, 2];
sim.drlt_cps = {[3], [1, 2, 3]};

shp.drlt_bnd = [1, 2];
shp.drlt_cps = {[1, 2, 3], [1, 2, 3]};

% Discretization parameters
deg = 3;
sub = 4;

shp.deg = [deg, deg];
shp.reg = [deg-1, deg-1];
shp.sub = [sub, sub];
shp.ngp = [deg+1, deg+1];

deg = 3;
sub = 1;

sim.deg = [deg, deg];
sim.reg = [deg-1, deg-1];
sim.sub = [sub, sub];
sim.ngp = [deg+1, deg+1];

% Geometry parameters
geo.rad = 10;
geo.len = 10;
geo.acv = 100*pi;

geo.nrb = half_cylinder_surface(geo.rad, geo.len);
shp.nrb = nrbrefine(geo.nrb, shp.deg, shp.reg, shp.sub);

init_dsn = shp.nrb.coefs(1:3, :, :);
init_dsn = init_dsn(:);

% Visualization parameters
vis.sub = [120, 40];

%% Initial design
shp.nrb = setup_nrb(init_dsn, shp.nrb);
sim.nrb = nrbrefine(shp.nrb, sim.deg, sim.reg, sim.sub);

sim = assemble_koiter(sim);
shp = assemble_shape(shp);

[obj, drv] = objective_koiter(sim, shp);
[ceq, deq] = constraint_volume(geo.acv, shp);

figure;
nrbplot(sim.nrb, vis.sub);
axis equal;
title('Initial design', 'FontSize', 12);

figure;
dfmplot(sim.sol, sim.nrb, sim.spc, vis.sub);
axis equal;
title('Deformed initial design', 'FontSize', 12);

fprintf('[Initial design]\n');
fprintf('Objective function value: %f\n', obj);
fprintf('Objective gradient norm: %f\n', shape_norm(shape_gradient(drv, shp), shp));
fprintf('Constraint function value: %f\n', ceq);
fprintf('Constraint gradient norm: %f\n', shape_norm(shape_gradient(deq, shp), shp));

%% Optimization
global glob_iter objc_vals;
glob_iter = -1;
objc_vals = [];

objective = @(dsn) problem_koiter_fmin_objective(dsn, shp, sim);
constraint = @(dsn) problem_koiter_fmin_constraint(dsn, geo, shp, sim);

options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');
options = optimoptions(options, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
options = optimoptions(options, 'MaxIterations', 25, 'MaxFunctionEvaluations', Inf, 'OutputFcn', @(dsn, optimValues, state) output_fmc(dsn, optimValues, state, shp));

figure;
tic;
[final_dsn, final_obj, ~, final_out] = fmincon(objective, init_dsn, [], [], [], [], [], [], constraint, options);
elapsed = toc;
fprintf('Time elapsed: %f\n', elapsed);

%% Final design
shp.nrb = setup_nrb(final_dsn, shp.nrb);
sim.nrb = nrbrefine(shp.nrb, sim.deg, sim.reg, sim.sub);

sim = assemble_koiter(sim);
shp = assemble_shape(shp);

[obj, drv] = objective_koiter(sim, shp);
[ceq, deq] = constraint_volume(geo.acv, shp);

output_fcn(final_dsn, obj, drv, drv, shp, vis);

%% Visualization
figure;
nrbplot(sim.nrb, vis.sub);
axis equal;
title('Final design', 'FontSize', 12);

figure;
dfmplot(sim.sol, sim.nrb, sim.spc, vis.sub);
axis equal;
title('Deformed final design', 'FontSize', 12);

fprintf('[Final design]\n');
fprintf('Objective function value: %f\n', obj);
fprintf('Objective gradient norm: %f\n', shape_norm(shape_gradient(drv, shp), shp));
fprintf('Constraint function value: %f\n', ceq);
fprintf('Constraint gradient norm: %f\n', shape_norm(shape_gradient(deq, shp), shp));

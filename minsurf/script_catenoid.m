% Problem parameters
shp.drlt_bnd = [1, 2, 3, 4];
shp.drlt_cps = {[2], [2], [1, 2, 3], [1, 2, 3]};

% Discretization parameters
deg = 3;
sub = 2;

shp.deg = [deg, deg];
shp.reg = [deg-1, deg-1];
shp.sub = [sub, sub];
shp.ngp = [deg+1, deg+1];

% Geometry parameters
geo.rad = 1;
geo.len = 1;

geo.nrb = cylinder_surface(geo.rad, geo.len);
shp.nrb = nrbrefine(geo.nrb, shp.deg, shp.reg, shp.sub);

init_dsn = shp.nrb.coefs(1:3, :, :);
init_dsn = init_dsn(:);

% Visualization parameters
vis.sub = [60, 10];

%% Initial design
shp.nrb = setup_nrb(init_dsn, shp.nrb);
shp = assemble_shape(shp);
[obj, drv] = objective_volume(shp);

figure;
nrbplot(shp.nrb, vis.sub);
axis equal;
title('Shape of initial design');

fprintf('[Initial design]\n');
fprintf('Objective function value: %f\n', obj);
fprintf('Objective gradient norm: %f\n', shape_norm(shape_gradient(drv, shp), shp));

%% Optimization
global glob_iter objc_vals;
glob_iter = 0;
objc_vals = [];

options.func_val_tol = 1e-6;
options.grad_nrm_tol = 1e-3;
options.min_step_size = 1e-6;
options.max_step_size = Inf;
options.init_step_size = 1e-2;
options.sigma = 1e-4;
options.tau = 0.5;
options.output_fcn = @(dsn, obj, dsc, nrm) output_fcn(dsn, obj, dsc, nrm, shp, vis);

problem = @(dsn) problem_minsurf(dsn, shp);

figure;
tic;
final_dsn = descent_method(problem, init_dsn, options);
elapsed = toc;
fprintf('Time elapsed: %f\n', elapsed);

%% Final design
shp.nrb = setup_nrb(final_dsn, shp.nrb);
shp = assemble_shape(shp);
[obj, drv] = objective_volume(shp);

%figure;
nrbplot(shp.nrb, vis.sub);
axis equal;
title('Shape of final design');

fprintf('[Final design]\n');
fprintf('Objective function value: %f\n', obj);
fprintf('Objective gradient norm: %f\n', shape_norm(shape_gradient(drv, shp), shp));

%% https://mathworld.wolfram.com/MinimalSurfaceofRevolution.html
% minimal objective value for geo.rad = 1 and geo.len = 1: 5.9917...
% catenary factor a = 0.8483

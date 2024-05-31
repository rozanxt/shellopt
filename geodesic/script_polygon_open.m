A = 1.0;
dt = 0.1;

n = 101;
c = [linspace(0, 1, n); zeros(1, n)];

[~, ~, cs, ~, ~, nc, ~, op] = compute_polygon_open(c, A);

t = linspace(-1, 1, 62);
a0 = zeros(1, n);
a0(21 : 80) = exp(-1 ./ (1 - t(2 : 61).^2));
an0 = repmat(a0, 2, 1) .* nc;
anvec0 = an0';
anvec0 = anvec0(:);
dof = [2 : (n - 1), (n + 2) : (2 * n - 1)];
ct = reshape(op(dof, dof) \ anvec0(dof), n - 2, 2)';
ct = [zeros(2, 1), ct, zeros(2, 1)];
a = a0;

figure;
for i = 1 : 300
    c = c + dt * ct;
    
    [dc, lc, cs, fc, ds, nc, kc, op] = compute_polygon_open(c, A);
    
    dsct = (ct(:, [2 : n, 1]) - ct(:, [n, 1 : (n - 1)])) ./ repmat(2 * cs, 2, 1);
    bv = (0.5 * kc) .* (A * dot(dsct, dsct) - dot(ct, ct)) - dot(dsct, ds) .* a;
    
    a = a + dt * bv;
    an = repmat(a, 2, 1) .* nc;
    anvec = an';
    anvec = anvec(:);
    dof = [2 : (n - 1), (n + 2) : (2 * n - 1)];
    ct = reshape(op(dof, dof) \ anvec(dof), n - 2, 2)';
    ct = [zeros(2, 1), ct, zeros(2, 1)];
    
    plot(c(1, :), c(2, :));
    hold on;
    quiver([c(1, :), c(1, 1)], [c(2, :), c(2, 1)], [ct(1, :), ct(1, 1)], [ct(2, :), ct(2, 1)]);
    hold off;
    axis equal;
    
    pause;
end

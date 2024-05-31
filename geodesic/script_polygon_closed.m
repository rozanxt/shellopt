A = 1.0;
dt = 0.1;
m = 100;
n = 200;

s = linspace(0, 2 * pi, n + 1);
c = [cos(s(1 : n)); sin(s(1 : n))];

[~, ~, cs, ~, ~, nc, ~, op] = compute_polygon_closed(c, A);

ct0 = [ones(1, n); zeros(1, n)];
ctvec0 = ct0';
ctvec0 = ctvec0(:);
a0 = dot(reshape(op * ctvec0, n, 2)', nc);
pn = a0 .* nc .* cs;

figure;
for i = 1 : m
    [~, ~, cs, ~, ~, nc, kc, op] = compute_polygon_closed(c, A);
    
    an = dot(nc, pn ./ cs) .* nc;
    anvec = an';
    anvec = anvec(:);
    ct = reshape(op \ anvec, n, 2)';
    c = c + dt * ct;
    
    dsct = (ct(:, [2 : n, 1]) - ct(:, [n, 1 : (n - 1)])) ./ (2 * cs);
    bn = (0.5 * kc) .* (A * dot(dsct, dsct) - dot(ct, ct)) .* nc .* cs;
    pn = pn + dt * bn;
    
    plot([c(1, :), c(1, 1)], [c(2, :), c(2, 1)]);
    hold on;
    quiver([c(1, :), c(1, 1)], [c(2, :), c(2, 1)], [ct(1, :), ct(1, 1)], [ct(2, :), ct(2, 1)]);
    hold off;
    axis equal;
    xlim([-1, 5]);
    ylim([-3, 3]);
    pause;
end

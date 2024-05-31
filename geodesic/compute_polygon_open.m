function [dc, lc, cs, fc, ds, nc, kc, op] = compute_polygon_open(c, A)
    n = size(c, 2);
    
    dc = c(:, [2 : n, 1]) - c(:, 1 : n);
    lc = sqrt(dc(1, :).^2 + dc(2, :).^2);
    cs = 0.5 * (lc([n, 1 : (n - 1)]) + lc(1 : n));
    fc = dc ./ repmat(lc, 2, 1);
    ds = (c(:, [2 : n, 1]) - c(:, [n, 1 : (n - 1)])) ./ repmat(lc([n, 1 : (n - 1)]) + lc(1 : n), 2, 1);
    
    nc = (fc(:, [n, 1 : (n - 1)]) + fc(:, 1 : n));
    nc = nc ./ repmat(sqrt(nc(1, :).^2 + nc(2, :).^2), 2, 1);
    nc = flip(nc);
    nc(1, :) = -nc(1, :);
    nc(:, 1) = [0; 1];
    nc(:, end) = [0; 1];
    
    kc = dot(fc(:, [n, 1 : (n - 1)]) - fc(:, 1 : n), nc);
    %kc = acos(dot(fc(:, [n, 1 : (n - 1)]), fc(:, 1 : n)));
    %kc = atan2(fc(2, 1 : n), fc(1, 1 : n)) - atan2(fc(2, [n, 1 : (n - 1)]), fc(1, [n, 1 : (n - 1)]));
    
    h1 = lc([n, 1 : (n - 1)])';
    h2 = lc(1 : n)';
    tp = [h1 .* (h1 + h2), -(h1 .* h2), h2 .* (h1 + h2)];
    tp = 2 ./ tp;
    lp = diag(tp(2 : end, 1), -1) + diag(tp(:, 2)) + diag(tp(1 : (end - 1), 3), 1);
    lp(1, end) = tp(1, 1);
    lp(end, 1) = tp(end, 3);
    op = eye(n) - A * lp;
    op = blkdiag(op, op);
end

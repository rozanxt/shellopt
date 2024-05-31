function nrb = unit_square
    nrb_pts = zeros(4, 2, 2);
    nrb_pts(:, 1, 1) = [0, 0, 1, 1]';
    nrb_pts(:, 2, 1) = [1, 0, 1, 1]';
    nrb_pts(:, 1, 2) = [0, 1, 1, 1]';
    nrb_pts(:, 2, 2) = [1, 1, 1, 1]';
    nrb_knt = {[0, 0, 1, 1], [0, 0, 1, 1]};
    nrb = nrbmak(nrb_pts, nrb_knt);
end

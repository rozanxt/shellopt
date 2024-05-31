function nrb = half_cylinder_surface(rad, len)
    nrb_pts = zeros(4, 4);
    nrb_pts(:, 1) = [rad, 0, 0, 1]';
    nrb_pts(:, 2) = [rad, rad, 0, 1]'/2;
    nrb_pts(:, 3) = [-rad, rad, 0, 1]'/2;
    nrb_pts(:, 4) = [-rad, 0, 0, 1]';
    nrb_knt = [0, 0, 0, 0.5, 1, 1, 1];
    nrb = nrbmak(nrb_pts, nrb_knt);
    nrb = nrbextrude(nrbtform(nrb, vecrotx(pi/2)), [0, len, 0]);
end

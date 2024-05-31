function nrb = nrbrefine(nrb, deg, reg, sub)
    nrb = nrbdegelev(nrb, max(deg-(nrb.order-1), 0));
    [~, ~, knt] = kntrefine(nrb.knots, sub-1, nrb.order-1, reg);
    nrb = nrbkntins(nrb, knt);
end

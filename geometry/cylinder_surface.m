function nrb = cylinder_surface(rad, len)
    nrb = nrbcirc(rad);
    nrb = nrbextrude(nrb, [0, 0, len]);
end

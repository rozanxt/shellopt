function [msh, nds, wts, rule] = setup_msh(nrb, ngp)
    rule = msh_gauss_nodes(ngp);
    [nds, wts] = msh_set_quad_nodes(nrb.knots, rule);
    msh = msh_cartesian(nrb.knots, nds, wts, geo_load(nrb));
end

function dfmplot(sol, nrb, spc, sub)
    cfg = geo_deform(sol, spc, geo_load(nrb));
    nrbplot(cfg.nurbs, sub);
end

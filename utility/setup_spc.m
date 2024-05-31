function spc = setup_spc(nrb, msh)
    spc = sp_nurbs(nrb, msh);
    spc = sp_vector(repmat({spc}, 1, msh.rdim), msh);
end

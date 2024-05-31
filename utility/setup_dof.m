function [free_dof, drlt_dof] = setup_dof(spc, drlt_bnd, drlt_cps)
    drlt_dof = [];
    for i = 1:numel(drlt_bnd)
        sd = drlt_bnd(i);
        cp = drlt_cps{i};
        for j = cp
            bd = spc.boundary(sd);
            drlt_dof = union(drlt_dof, bd.dofs(bd.comp_dofs{j}));
        end
    end
    free_dof = setdiff(1:spc.ndof, drlt_dof);
end

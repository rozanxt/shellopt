function coefs = greville(knots, order)
    coefs = cell(1, length(knots));
    for i = 1:length(knots)
        knot = knots{i};
        ordi = order(i);
        dofs = length(knot)-ordi;
        coef = zeros(1, dofs);
        for j = 1:ordi-1
            coef = coef+knot(1+j:dofs+j);
        end
        coefs{i} = coef;
    end
end

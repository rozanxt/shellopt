function design = penalty_method(problem, design, options, penalty)
    for pen = penalty
        design = descent_method(@(dsn) problem(dsn, pen), design, options);
    end
end

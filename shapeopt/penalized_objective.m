function [obj, drv] = penalized_objective(obj, drv, ceq, deq, pen)
    obj = obj+pen*ceq^2;
    drv = drv+2*pen*ceq*deq;
end

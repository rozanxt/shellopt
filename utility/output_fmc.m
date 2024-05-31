function stop = output_fmc(dsn, optimValues, state, shp)
    obj = optimValues.fval;
    global glob_iter objc_vals;
    objc_vals = [objc_vals, obj];
    nrbplot(setup_nrb(dsn, shp.nrb), [120, 40]);
    axis equal;
    %xlim([-10.5, 10.5]);
    %ylim([-0.5, 10.5]);
    %zlim([0, 11.5]);
    title(sprintf('Iteration: %d, Objective: %f', glob_iter, obj), 'FontSize', 12);
    %exportgraphics(gcf, ['iteration-', num2str(glob_iter), '.png'], 'Resolution', 200);
    stop = false;
    glob_iter = glob_iter + 1;
    pause(0.001);
end

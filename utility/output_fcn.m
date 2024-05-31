function output_fcn(dsn, obj, ~, nrm, shp, vis)
    global glob_iter objc_vals;
    objc_vals = [objc_vals, obj];
    nrbplot(setup_nrb(dsn, shp.nrb), vis.sub);
    axis equal;
    %xlim([-10.5, 10.5]);
    %ylim([-0.5, 10.5]);
    %zlim([0, 11.5]);
    title(sprintf('Iteration: %d, Objective: %f', glob_iter, obj), 'FontSize', 12);
    %exportgraphics(gcf, ['iteration-', num2str(glob_iter), '.png'], 'Resolution', 200);
    glob_iter = glob_iter + 1;
    pause(0.001);
end

function prev_dsn = descent_method(problem, init_dsn, options)
    prev_dsn = init_dsn;
    step_size = options.init_step_size;
    
    num_iter = 0;
    fcn_eval = 0;
    
    terminate = false;
    while ~terminate
        [prev_obj, desc_dir, grad_nrm] = problem(prev_dsn);
        fcn_eval = fcn_eval+1;
        
        options.output_fcn(prev_dsn, prev_obj, desc_dir, grad_nrm);
        
        if grad_nrm < options.grad_nrm_tol
            fprintf('Stationary point found (gradient norm criterion).\n');
            break;
        end
        
        next_dsn = prev_dsn+step_size*desc_dir;
        next_obj = problem(next_dsn);
        fcn_eval = fcn_eval+1;
        
        if next_obj < prev_obj-step_size*options.sigma*grad_nrm^2
            while true
                step_size = step_size/options.tau;
                
                if step_size > options.max_step_size
                    break;
                end
                
                temp_dsn = prev_dsn+step_size*desc_dir;
                temp_obj = problem(temp_dsn);
                fcn_eval = fcn_eval+1;
                
                if temp_obj < prev_obj-step_size*options.sigma*grad_nrm^2
                    next_dsn = temp_dsn;
                    next_obj = temp_obj;
                else
                    break;
                end
            end
        end
        
        step_size = step_size*options.tau;
        
        while next_obj > prev_obj-step_size*options.sigma*grad_nrm^2
            step_size = step_size*options.tau;
            
            if step_size < options.min_step_size
                fprintf('No decrease along descent direction (step size criterion).\n');
                next_dsn = prev_dsn;
                terminate = true;
                break;
            end
            
            next_dsn = prev_dsn+step_size*desc_dir;
            next_obj = problem(next_dsn);
            fcn_eval = fcn_eval+1;
        end
        
        if abs(next_obj-prev_obj) < options.func_val_tol*(1+abs(prev_obj))
            fprintf('Stationary point found (function value criterion).\n');
            terminate = true;
        end
        
        prev_dsn = next_dsn;
        num_iter = num_iter+1;
    end
    
    fprintf('Number of iterations: %d, Number of function evaluations: %d\n', num_iter, fcn_eval);
end

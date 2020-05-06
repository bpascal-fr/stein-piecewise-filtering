function [xfilter,lambda, crit]= bfgs_sugar(data, param, op, dopl, prox, dprox, objective)
    % Automated tuning of hyperparamters
    % Performs a quasi-Newton descent minimizing SURE
    
    
    % Evaluate the standard deviation of the noise
    if strcmp(param.type,'1D')
        Nl = size(data,1);
        if isfield(param,'sigma')
            sigma = param.sigma.*ones(1,Nl);
        else
            sigma = zeros(1,Nl);
            for nl = 1:Nl
                [~,cH] = dwt(data(nl,:),'db1');
                C = abs(cH);
                sigma(nl) = median(C)/0.6745;
            end
        end
    elseif strcmp(param.type,'2D')
        if isfield(param,'sigma')
            sigma = param.sigma;
        else
            [~,cH,cV,cD] = dwt2(data,'db1');
            C = abs([cH(:) ; cV(:); cD(:)]);
            sigma = median(C)/0.6745;
        end
    else
        error('Only (multivariate) signals of type 1D or images of type 2D allowed')
    end
    
    % Compute Finite Difference Monte Carlo parameters
    % sure.eps: finite difference step
    % sure.delta: Monte Carlo unitary Gaussian vector
    if isfield(param,'fdmc')
        fdmc = param.fdmc;
    else
        if isfield(param,'nu')
            fdmc.eps = param.nu;
        else
            fdmc.eps = 2*max(sigma(:))/numel(data)^.3;
        end
        fdmc.delta = randn(size(data));
        fdmc.sigma = sigma;
    end
    
    % Initialize the quasi-Newton algorithm
    
    % Initial hyperparameter
    if param.lambda == - 1
        Nl = 1;
    else
        if strcmp(param.type,'1D')
            Nl = size(data,1);
        else
            error('Multivariate only for signals')
        end
    end
    lambda_in = zeros(Nl,1);
    if strcmp(param.type,'1D')
        for nl = 1:Nl
            lambda_in(nl) = numel(data)*sigma(nl)^2/(4*objective.regularization(op.direct(data(nl,:),1),1));
        end
    else
        lambda_in = numel(data)*sigma^2/(4*objective.regularization(op.direct(data,1),1));
    end
    opts.x0 = lambda_in;
    % Initial Hessian
    set_init(struct);
    set_crits(struct);
    [~,sugar] = sure_sugar(data, param, op, dopl, prox, dprox, objective, fdmc,lambda_in);
    if isfield(param,'alpha')
        alpha = param.alpha;
    else
        alpha = 0.9;
    end
    opts.H0 = diag(abs(alpha*lambda_in)./abs(sugar));
    
    
    % Define SURE and SUGAR as a function
    sure_sugar_fun = @(lambda) sure_sugar(data, param, op, dopl, prox, dprox, objective, fdmc, lambda);
    
    % Run BFGS GRANSO algorithm
    opts.print_level = 1;
    opts.quadprog_info_msg = false;
    opts.prescaling_info_msg = false;
    opts.maxit = 20;
    soln    = granso(Nl,sure_sugar_fun,@positivityConstraint,[],opts);
    lambda = soln.final.x;
    
    init_PD = get_init;
    xfilter = init_PD.x;
    crit_PD = get_crits;
    crit_PD.lambda = cell2mat(crit_PD.lambda);
    crit = crit_PD;
    
    % Define positivity constraints
    function [ci,ci_grad] = positivityConstraint(lambda)
        % Impose that lambda >= 1e-2 * lambda_in (so that lambda > 0)
        ci      = 1e-2 * lambda_in - lambda;
        ci_grad = -eye(Nl);
    end
end
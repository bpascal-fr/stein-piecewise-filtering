function [xfilter,lambda,crit,gap] = nonlinear_filtering(data,type,filter_def,lambda,norm_data,norm_reg,computation,param)
    %
    % Perform linear filtering in 1D
    %
    % [xfilter,lambda] = nonlinear_filtering(data,filter_def,lambda,norm_data,norm_reg,computation)
    % data         : data to be filtered
    % filter_def   : filter can either be defined as 'gradient' or 'laplacian'
    %               or on a form of the PSF (i.e. [1/2, -1/2]) or the OTF
    % lambda       : regularisation parameter (the larger, the smoother)
    %               default value N^(1/2)*sigma/4
    % norm_data    : norm of the data-term can be 1 (impulsional noise) or 2
    % xfilter      : denoised signal
    % crit         : objective function value
    % param.sigma  : knwon standard deviation of the noise (optional)
    % param.tol    : tolerance on objective increments to stop the minimization (optional)
    % param.nu     : finite difference step for differentiation w.r.t. observations (optional)
    % param.alpha  : parameter for initializing the Hessian (optional, default 0.9)
    % param.fdmc   : parameters for Finite Difference Monte Carlo SURE
    %
    % minimization problem:
    % --------------------
    %       min_x   1/2|| data - x ||_2^2 + lambda ||h*x||_1
    %
    % exemple:
    % --------
    % >>  [xfilter,lambda] = nonlinear_filtering_1D(data,[-1/2,1/2],lambda)
    %
    % Implementation N. PUSTELNIK, ENS Lyon
    % June 2019
    
    if ~exist('norm_data') || isempty(norm_data) %#ok<EXIST>
        norm_data = 2;
    end
    
    if ~exist('norm_reg') || isempty(norm_reg)%#ok<EXIST>
        norm_reg = 1;
    end
    
    
    if strcmp('type','1D')
        [data_n,data_m] = size(data);
        if data_n>data_m
            data = data';
        end
        clear data_n data_m;
    end
    
    if ~exist('filter_def') || isempty(filter_def) %#ok<EXIST>
        filter_def      = 'gradient'; % Type of filter in the regularization
    end
    
    if ~exist('computation')|| isempty(computation) %#ok<EXIST>
        computation      = 'fourier'; % Type of filter in the regularization
    end
    
    if nargin < 8
        param = struct;
    end
    
    %% Objective function design
    param.type      = type ;         % or '2D' (multicomponent 1D is included in '1D')
    % Iteration number
    if ~isfield(param,'iter')
        if strcmp(param.type,'1D')
            param.iter      = 1e7;
        else
            param.iter      = 1e5;
        end
    end
    if ~isfield(param,'tol')
        if strcmp(param.type,'1D')
            param.tol       = 10^(-3);
        else
            param.tol       = 10^(-5);
        end
    end
    param.mu        = 1;
    param.lambda    = lambda;        % Regularization parameter
    param.p         = norm_data;
    param.q         = norm_reg;
    if param.p==1
        prox.fidelity = @(y,data,tau) prox_L1(y-data,tau)+data;
        objective.fidelity =  @(y,data) sum(abs(y(:)-data(:)));
    elseif param.p==2
        param.mu=1;
        prox.fidelity = @(y,data,tau) prox_L2(y-data,tau)+data;
        objective.fidelity =  @(y,data) 1/2*sum(abs(y(:)-data(:)).^2);
    end
    if param.q==1
        prox.regularization = @(y,tau) prox_L1(y,tau);
        objective.regularization =  @(y,tau) tau*sum(abs(y(:)));
    elseif param.q==2
        prox.regularization = @(y,tau) prox_L2(y,tau);
        objective.regularization =  @(y,tau) tau/2*sum(abs(y(:)).^2);
    elseif param.q==12
        if strcmp(param.type,'1D')
            prox.regularization = @(y,tau) prox_L12(y,tau);
            objective.regularization =  @(y,tau) tau*sum(sqrt(sum(y.^2,1)));
        else
            prox.regularization = @(y,tau) prox_L12_2D(y,tau);
            objective.regularization =  @(y,tau) tau*sum(sum(sqrt(sum(y.^2,3))));
        end
    end
    
    
    
    
    if numel(lambda) == 1
        if lambda == -1 || lambda == -2
            % Primal-dual with regularization parameter estimated with SUGAR
            op.direct         = @(x,lambda)opL_auto(x,filter_def,computation,param,lambda);
            op.normL          = @(x,lambda)normL2_auto(x,filter_def, computation,param,op,lambda);
            dopl.direct       = @(x,lambda,nl)dopL(x,filter_def,computation,param,lambda,nl);
            op.adjoint        = @(x,lambda)opLadj_auto(x,filter_def,computation,param,lambda);
            dopl.adjoint      = @(x,lambda,nl)dopLadj(x,filter_def,computation,param,lambda,nl);
            if param.p==1
                dprox.fidelity = @(y,data,dy,tau) dprox_L1(y-data,dy,tau);
            elseif param.p==2
                param.mu=1;
                dprox.fidelity = @(y,data,dy,tau) prox_L2(dy,tau);
            end
            if param.q==1
                dprox.regularization = @(y,dy,tau) dprox_L1(y,dy,tau);
            elseif param.q==2
                dprox.regularization = @(y,dy,tau) prox_L2(dy,tau);
            elseif param.q==12
                if strcmp(param.type,'1D')
                    dprox.regularization = @(y,dy,tau) dprox_L12(y,dy,tau);
                else
                    dprox.regularization = @(y,dy,tau) dprox_L12_2D(y,dy,tau);
                end
            end
            
            if isfield(param,'sigma')
                if ~(numel(param.sigma) == 1 || numel(param.sigma) == size(data,1))
                    error('param.sigma should be a scalar or a vector of length equal to the number of components')
                end
            end
            [xfilter,lambda,crit]  = bfgs_sugar(data, param, op, dopl, prox, dprox, objective);
        else
            % Primal-dual with regularization parameter fixed by user
            op.direct       = @(x)opL(x,filter_def,computation,param);
            op.adjoint      = @(x)opLadj(x,filter_def,computation,param);
            op.normL        = @(x)normL2(x,filter_def, computation,param,op);
            param.normL     = op.normL(data);
            [xfilter,obj,gap]  = PD_ChambollePock(data, param, op, prox, objective);
            crit.objective = obj;
            crit.gap = gap;
        end
    else
        % Primal-dual with regularization parameter fixed by user
        op.direct       = @(x)opL(x,filter_def,computation,param);
        op.adjoint      = @(x)opLadj(x,filter_def,computation,param);
        op.normL        = @(x)normL2(x,filter_def, computation,param,op);
        param.normL     = op.normL(data);
        [xfilter,obj,gap]  = PD_ChambollePock(data, param, op, prox, objective);
        crit.objective = obj;
        crit.gap = gap;
    end
    
end

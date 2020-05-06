function [param,op,dopl,prox,dprox,objective,fdmc] = nonlinear_filtering_prep(data,type,filter_def, lambda, norm_data,norm_reg,computation,param)
    %
    % Perform linear filtering in 1D
    %
    % [xfilter,lambda] = nonlinear_filtering(data,filter_def,lambda,norm_data,norm_reg,computation)
    % data        : data to be filtered
    % filter_def  : filter can either be defined as 'gradient' or 'laplacian'
    %               or on a form of the PSF (i.e. [1/2, -1/2]) or the OTF
    % lambda      : regularisation parameter type -1 (global) -2 (one per
    %               component, only for signals)
    % norm_data   : norm of the data-term can be 1 (impulsional noise) or 2
    % xfilter     : denoised signal
    % crit        : objective function value
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
    if lambda == -1
        param.lambda = 1;
    else
        if strcmp(type,'1D')
            param.lambda = ones(1,size(data,1));
        else
            error('Multiple lambda only available for signals i.e. type = ''1D''')
        end
    end
        
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
    
    param.sigma = sigma;
    
    fdmc.eps = 2*max(param.sigma(:))/numel(data)^.3;
    fdmc.delta = randn(size(data));
    fdmc.sigma = param.sigma;
    
    
end

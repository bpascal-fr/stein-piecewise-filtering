function [x_nlf_grid,lambda_grid,sure_vals,lambda_vals,x_grid,lambda_ref,fdmc] = grid_sure(data,type,filter_def,lambda,norm_reg,computation,param,grid,NL)
    
    % Compute Stein Unbiased Risk Estimate on a grid of hyperparameters
    % Inputs:  - data: signal or image to be analyzed
    %          - type: '1D' or '2D'
    %          - filter_def: filter can either be defined as 'gradient' or 'laplacian'
    %            or on a form of the PSF (i.e. [1/2, -1/2]) or the OTF
    %          - lambda: regularisation parameter type -1 (global) -2 (one per
    %            component, only for signals)
    %          - param: prior knowledge (sigma, lambda_ref)
    %          - grid: type of grid 'lin' (linearly spaced) or 'log'
    %            (logarithmically spaced)
    %          - NL: number of explored values of hyperparameter
    % Outputs: - x_nlf_grid: best estimate found in the grid (i.e. lowest SURE)
    %          - lambda_grid: corresponding best hyperparameter
    %          - sure_vals: SURE at each point of the grid
    %          - lambda_vals: values of hyperparameter explored
    %          - x_grid: estimate at each point of the grid
    %          - lambda_ref: center of the grid
    %          - fdmc: parameters for Finite Difference Monte Carlo
    % Comment: SURE is not defined for norm_data = 1 /!\, norm_data = 2
    % (mandataroy)
    %   
    % Implementation B. Pascal, ENS Lyon
    % June 2019
    
    
    [param,op,dopl,prox,dprox,objective,fdmc] = nonlinear_filtering_prep(data,type,filter_def,lambda,2,norm_reg,computation,param);
    Nl = size(data,1);
    
    if isfield(param,'sigma')
        sigma = param.sigma;
    else
        if Nl ==1
            sigma = std(data(:));
        else
            if Nl == 2
                sigma = zeros(1,Nl);
                for nl = 1:Nl
                    sigma(nl) = std(data(nl,:));
                end
            else
                error('Grid search only for 1 or 2 hyperparameters. More is too costly')
            end
        end
    end
    
    if isfield(param,'lambda_ref')
        lambda_ref = param.lambda_ref;
    else
        if lambda == -1
            lambda_ref = numel(data)*sigma^2/(4*objective.regularization(op.direct(data,1),1));
        else
            if Nl == 2
                lambda_ref = zeros(1,Nl);
                for nl = 1:Nl
                    lambda_ref(nl) = numel(data(nl,:))*sigma(nl)^2/(4*objective.regularization(op.direct(data(nl,:),1),1));
                end
            else
                error('Grid search only for 1 or 2 hyperparameters. More is too costly')
            end
        end
    end
    
    
    if strcmp(grid,'lin')
        lambda_ref = 5*lambda_ref/NL;
        LBD = lambda_ref'*(1:NL);
    elseif strcmp(grid,'log')
        LBD = lambda_ref'*logspace(-1,3,NL);
    else
        error('Wrong type of grid, choose ''lin'' or ''log''')
    end
   
    
    if numel(lambda_ref) == 1
        sure_vals = zeros(size(LBD));
        lambda_vals = LBD;
        x_grid = cell(1,NL);
        set_init(struct);
        set_crits(struct);
        for iL = 1:NL
            [sure_vals(iL), ~, x_grid{iL}] = sure_sugar(data, param, op, dopl, prox, dprox, objective, fdmc, LBD(iL));
            clc
            disp(['Grid search ',num2str(100*iL/NL,3),'% achieved'])
        end
        
        lambda_grid = LBD(sure_vals == min(sure_vals(:)));
        x_nlf_grid = x_grid{sure_vals == min(sure_vals(:))};
    else
        sure_vals = zeros(NL,NL);
        lambda_vals = LBD;
        x_grid = cell(NL,NL);
        set_init(struct);
        set_crits(struct);
        for iL_1 = 1:NL
            for iL_2 = 1:NL
                [sure_vals(iL_1,iL_2), ~, x_grid{iL_1,iL_2}] = sure_sugar(data, param, op, dopl, prox, dprox, objective, fdmc, [LBD(1,iL_1),LBD(2,iL_2)]);
                clc
                disp(['Grid search ',num2str(100*((iL_1-1)*NL+iL_2)/NL^2,3),'% achieved'])
            end
        end
        
        [ind_1,ind_2] = find(sure_vals == min(sure_vals(:)));
        lambda_grid = [LBD(1,ind_1), LBD(2,ind_2)];
        x_nlf_grid = x_grid{ind_1,ind_2};
    end
end
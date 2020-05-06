function [sure,sugar,x] = sure_sugar(data, param, op, dopl, prox, dprox, objective, fdmc, lambda)
% Compute SURE and SUGAR from iteratively differentiated algorithm dPD_ChambollePock
% for fixed regularization parameter lambda
% 
% Implementation B. PASCAL, ENS Lyon
% from ?
% - Deledalle, C.-A. and Vaiter, S. and Fadili, J. and Peyré, G.: Stein
% Unbiased GrAdient estimator of the Risk (SUGAR) for multiple parameter 
% selection. SIAM J. Imaging Sci. (2014)
% and 
% -  Pascal?B. and Vaiter S. and Pustelnik N. and Abry P.: Automated 
% data-driven selection of the hyperparameters for Total-Variation based 
% texture segmentation. Preprint arXiv:2004.09434 (2020)
% April 2020


    param.lambda  = lambda';
    
    % Run differentiated primal-dual
    [x,dx,Ex, Edx,obj,gap]=dPD_ChambollePock(data, param, op, dopl, prox, dprox, objective, fdmc);
    
    % Compute SURE
    if numel(fdmc.sigma) == 1
        sure = norm(x-data,'fro')^2 + 2*fdmc.sigma^2*sum((Ex(:)-x(:)).*fdmc.delta(:))/fdmc.eps - fdmc.sigma^2*numel(data);
    else
        Nl = size(data,1);
        sure = norm(x-data,'fro')^2;
        for nl = 1:Nl
            sure = sure + 2*fdmc.sigma(nl)^2*sum((Ex(nl,:)-x(nl,:)).*fdmc.delta(nl,:))/fdmc.eps - fdmc.sigma(nl)^2*numel(data(nl,:));
        end
    end
    
    % Compute SUGAR
    sugar = zeros(numel(param.lambda),1);
    for nlbd = 1:numel(param.lambda)
        if numel(fdmc.sigma) ==1
            sugar(nlbd) = 2*sum(dx{nlbd}(:).*(x(:)-data(:))) + 2*fdmc.sigma^2*sum((Edx{nlbd}(:)-dx{nlbd}(:)).*fdmc.delta(:))/fdmc.eps;
        else
            sugar(nlbd) = 0;
            for nl = 1:Nl
                sugar(nlbd) = sugar(nlbd) + 2*sum(dx{nlbd}(nl,:).*(x(nl,:)-data(nl,:))) + 2*fdmc.sigma(nl)^2*sum((Edx{nlbd}(nl,:)-dx{nlbd}(nl,:)).*fdmc.delta(nl,:))/fdmc.eps;
            end
        end
    end
    
    crit_PD = get_crits;
    if isfield(crit_PD,'neval')
        crit_PD.neval = crit_PD.neval + 1;
    else
        crit_PD.neval = 1;
    end
    crit_PD.objective{crit_PD.neval} = obj;
    crit_PD.gap{crit_PD.neval} = gap;
    crit_PD.lambda{crit_PD.neval} = lambda;
    crit_PD.sure(crit_PD.neval) = sure;
    set_crits(crit_PD);
    
end
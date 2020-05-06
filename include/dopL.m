function dxt = dopL(x,filter_def, computation, param, lambda, nl)
% Derivative of the linear operator associated with the filter in the prior 
% with respect to the nl-th regularization parameter 
% 
% 
% Implementation B. PASCAL, ENS Lyon
% April 2020

    dlambda = zeros(size(lambda));
    dlambda(nl) = 1;
    param.lambda = dlambda;
    dxt = opL(x,filter_def, computation,param);

end

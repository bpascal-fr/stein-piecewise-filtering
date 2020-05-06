
function [x,obj,gap]=PD_ChambollePock(data, param, op, prox, objective)
% Primal-dual algorithm by Chambolle and Pock handling strong convexity 
% when possible. 
% Chambolle, A., Pock, T.: A first-order primal-dual algorithm for convex 
% problems with applications to imaging. J. Math. Imag. Vis. 40(1), 
% 120?145 (2011)
% 
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% April 2020

    
    %% Fixing Proximal Parameters
    gamma = 0.99;
    tau = gamma/sqrt(param.normL);
    sig = gamma/sqrt(param.normL);
    if tau*sig*param.normL>1
        disp('ERROR');
    end
    theta=1;
    
    %% Initializing variables
    
    x = zeros(size(data));
    y = op.direct(x);
    x0 = x;
    bx = x;
    
    %% Criterion of convergence
    
    obj = zeros(1,param.iter);
    gap = obj;
    i = 0;
    gapc = param.tol + 1;
    
    %% Algorithm
    while (gapc > param.tol)&&(i < param.iter)
        i = i+1;
        
        %for i=1:param.iter
        %Update of primal variable
        tmp = y + sig*op.direct(bx);
        y = tmp - sig*prox.regularization(tmp/sig, 1/sig);
        
        %Update of dual variable
        x = prox.fidelity(x0 - tau * op.adjoint(y),data,tau);
        
        %Update of the descent steps
        if param.mu>=1
            theta = (1+2*param.mu*tau)^(-1/2);
            tau = theta*tau;
            sig=sig/theta;
        end
        
        %Update dual auxiliary variable
        bx = x + theta*(x - x0);
        x0 = x;
        obj(i) = objective.fidelity(x,data) + objective.regularization(op.direct(x),1);
        P = prox.regularization(y,1);
        dual =  objective.fidelity(op.adjoint(y),data) - 1/2*norm(data,'fro')^2 + norm(P(:),'fro');
        gap(i) = (obj(i) + dual)/(abs(obj(i)) + abs(dual));
        gapc = gap(i);
        
    end
    
    obj = obj(1:i);
    gap = gap(1:i);
end


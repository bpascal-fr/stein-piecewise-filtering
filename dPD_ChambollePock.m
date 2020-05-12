%**************************************************************************
% authors: Barbara Pascal                                                 *
% institution: Laboratoire de Physique de l'ENS de Lyon                   *
% date: March 17 2020                                                     *
% License CeCILL-B                                                        *
%**************************************************************************

function [x,dx,Ex,Edx,obj,gap]=dPD_ChambollePock(data, param, op, dopl, prox, dprox, objective, fdmc)
    
    
    Edata  = data + fdmc.eps*fdmc.delta;
    
    %% Fixing Proximal Parameters
    gamma = 0.99;
    lambda = param.lambda;
    normL = op.normL(data,lambda);
    tau = gamma/sqrt(normL);
    sig = gamma/sqrt(normL);
    if tau*sig*normL>1
        disp('ERROR');
    end
    theta=1;
    %% Initializing variables
    Nl = numel(param.lambda);
    init_PD = get_init;
    
    if isfield(init_PD,'x')
        dx = init_PD.dx;
        x = init_PD.x;
        Edx = init_PD.Edx;
        Ex = init_PD.Ex;
        dy = init_PD.dy;
        y = init_PD.y;
        Edy = init_PD.Edy;
        Ey = init_PD.Ey;
        dx0 = init_PD.dx0;
        x0 = init_PD.x0;
        Edx0 = init_PD.Edx0;
        Ex0 = init_PD.Ex0;
        dbx = init_PD.dbx;
        bx = init_PD.bx;
        Edbx = init_PD.Edbx;
        Ebx = init_PD.Ebx;
    else
        dx = cell(1,Nl);
        for nl = 1:Nl
            dx{nl} = zeros(size(data));
        end
        x = zeros(size(data));
        Edx = cell(1,Nl);
        for nl = 1:Nl
            Edx{nl} = zeros(size(data));
        end
        Ex = zeros(size(data));
        
        dy = cell(1,Nl);
        for nl = 1:Nl
            dy{nl} =  op.direct(dx{nl},lambda) + dopl.direct(x,lambda,nl);
        end
        y = op.direct(x,lambda);
        Edy = cell(1,Nl);
        for nl = 1:Nl
            Edy{nl} =  op.direct(Edx{nl},lambda) + dopl.direct(Ex,lambda,nl);
        end
        Ey = op.direct(Ex,lambda);
        
        dx0 = dx;
        x0 = x;
        Edx0 = Edx;
        Ex0 = Ex;
        
        dbx = dx;
        bx = x;
        Edbx = Edx;
        Ebx = Ex;
    end
    dtmp = cell(1,Nl);
    Edtmp = cell(1,Nl);
    
    %% Criterion of convergence
    
    obj = zeros(1,param.iter);
    dual = obj;
    gap = obj;
    i = 0;
    gapc = param.tol + 1;
    
    
    %% Algorithm
    while (gapc > param.tol)&&(i < param.iter)
        i = i+1;
        %for i=1:param.iter
        %Update of primal variable
        for nl = 1:Nl
            dtmp{nl} = dy{nl} + sig*op.direct(dbx{nl},lambda) + sig*dopl.direct(bx,lambda,nl);
        end
        tmp = y + sig*op.direct(bx,lambda);
        for nl = 1:Nl
            Edtmp{nl} = Edy{nl} + sig*op.direct(Edbx{nl},lambda) + sig*dopl.direct(Ebx,lambda,nl);
        end
        Etmp = Ey + sig*op.direct(Ebx,lambda);
        
        for nl = 1:Nl
            dy{nl} = dtmp{nl} - dprox.regularization(tmp/sig, dtmp{nl}, 1/sig);
        end
        y = tmp - sig*prox.regularization(tmp/sig, 1/sig);
        for nl = 1:Nl
            Edy{nl} = Edtmp{nl} - dprox.regularization(Etmp/sig, Edtmp{nl}, 1/sig);
        end
        Ey = Etmp - sig*prox.regularization(Etmp/sig, 1/sig);
        
        %Update of dual variable
        for nl = 1:Nl
            dx{nl} = dprox.fidelity(x0, data, dx0{nl} - tau * op.adjoint(dy{nl},lambda) - tau*dopl.adjoint(y,lambda,nl),tau);
        end
        x  = prox.fidelity(x0 - tau * op.adjoint(y,lambda),data,tau);
        for nl = 1:Nl
            Edx{nl} = dprox.fidelity(Ex0, Edata, Edx0{nl} - tau * op.adjoint(Edy{nl},lambda) - tau*dopl.adjoint(Ey,lambda,nl),tau);
        end
        Ex = prox.fidelity(Ex0 - tau * op.adjoint(Ey,lambda),Edata,tau);
        
        %Update of the descent steps
        if param.mu>=0
            theta = (1+2*param.mu*tau)^(-1/2);
            tau = theta*tau;
            sig=sig/theta;
        end
        
        %Update dual auxiliary variable
        for nl = 1:Nl
            dbx{nl} = dx{nl} + theta*(dx{nl} - dx0{nl});
        end
        bx = x + theta*(x - x0);
        for nl = 1:Nl
            Edbx{nl} = Edx{nl} + theta*(Edx{nl} - Edx0{nl});
        end
        Ebx = Ex + theta*(Ex - Ex0);
        
        for nl = 1:Nl
            dx0{nl} = dx{nl};
        end
        x0 = x;
        for nl = 1:Nl
            Edx0{nl} = Edx{nl};
        end
        Ex0 = Ex;
        
        
        obj(i) = objective.fidelity(x,data) + objective.regularization(op.direct(x,lambda),1);
        P = prox.regularization(y,1);
        dual(i) =  objective.fidelity(op.adjoint(y,lambda),data) - 1/2*norm(data,'fro')^2 + norm(P(:),'fro');
        gap(i) = (obj(i) + dual(i))/(abs(obj(i)) + abs(dual(i)));
        gapc = gap(i);
        
    end
    
    obj = obj(1:i);
    dual = dual(1:i);
    gap = gap(1:i);


    
    init_PD.dx = dx;
    init_PD.x = x;
    init_PD.Edx = Edx;
    init_PD.Ex = Ex;
    init_PD.dy = dy;
    init_PD.y = y;
    init_PD.Edy = Edy;
    init_PD.Ey = Ey;
    init_PD.dx0 = dx0;
    init_PD.x0 = x0;
    init_PD.Edx0 = Edx0;
    init_PD.Ex0 = Ex0;
    init_PD.dbx = dbx;
    init_PD.bx = bx;
    init_PD.Edbx = Edbx;
    init_PD.Ebx = Ebx;
    init_PD.crit = obj;
    init_PD.dual = dual;
    init_PD.gap = gap;
    
    set_init(init_PD);
end



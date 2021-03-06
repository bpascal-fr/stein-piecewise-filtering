function [dProx] = dprox_L12(P,dP,gamma)
% Derivative of the proximity operator of the l12 norm with respect to its
% first argument P, applied to dP
% dProx = (d/dP) prox_{gamma || .||_{1,2}}(P) [dP]
% 
% Implementation B. PASCAL, ENS Lyon
% April 2020


[n,~] = size(P);
if n==1
    dProx = dprox_L1(P,dP,gamma);
else
    nP=(sum(P.^2,1)).^(1/2);
    dProx = zeros(size(P));
    Proj = sum(dP(:,nP >gamma).*P(:,nP >gamma),1)./nP(nP >gamma).^2;
    ProjP = dP(:,nP >gamma) - Proj.*P(:,nP >gamma);
    dProx(:,nP >gamma) = dP(:,nP > gamma) - gamma*ProjP./nP(nP > gamma);
end

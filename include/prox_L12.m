function [Prox] = prox_L12(P,gamma)
% Proximity operator of the l12 norm 
% Prox = prox_{gamma || .||_{1,2}}(P)
% 
%
% Implementation B. PASCAL, ENS Lyon
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% April 2020

[n,~] = size(P);
if n==1
    Prox = prox_L1(P,gamma);
else
    temp=(sum(P.^2,1)).^(1/2);
    ind=find(temp>gamma);
    Prox = zeros(size(P));
    Prox(:,ind)=kron((1-gamma./temp(ind)),ones(n,1)).*P(:,ind);
end

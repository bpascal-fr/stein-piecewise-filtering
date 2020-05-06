function dProx = dprox_L12_2D(P,dP,gamma)
% Derivative of the proximity operator of the l12 norm with respect to its
% first argument P, applied to dP
% dProx = (d/dP) prox_{gamma || .||_{1,2}}(P) [dP]
%
% Implementation B. PASCAL, ENS Lyon
% April 2020

    [n1,n2,n] = size(P);
    dP = reshape(dP,n1*n2,n);
    P = reshape(P,n1*n2,n);

    dProx = dprox_L12(P,dP,gamma);
    dProx = reshape(dProx,n1,n2,n);
end

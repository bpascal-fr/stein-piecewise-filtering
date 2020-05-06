function Prox = prox_L12_2D(P,gamma)
% Proximity operator of the l12 norm 
% Prox = prox_{gamma || .||_{1,2}}(P)
% 
%
% Implementation B. PASCAL, ENS Lyon
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% April 2020

    [n1,n2,n] = size(P);
    P = reshape(P,n1*n2,n)';

    Prox = prox_L12(P,gamma);
    Prox = Prox';
    Prox = reshape(Prox,n1,n2,n);
end

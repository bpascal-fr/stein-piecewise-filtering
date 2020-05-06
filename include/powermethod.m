function beta=powermethod(x,op)

% Power method to estimate ||op||^2
%
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% June 2019

rhon=1+1e-6;
rhon1(1)=0;
[n,m] = size(x);
xn  = randn(n,m);
xn1 = xn;
k=1;
while abs(rhon1(k)-rhon)/rhon1(k) >= 1e-4
   xn  = xn1; 
   xn1 = op.adjoint(op.direct(xn));

   rhon=rhon1(k);
   k=k+1;
   rhon1(k) = sum(sum(sum(xn1.^2)))/sum(sum(sum(xn.^2)));
   fprintf('Norm^2 = %3.8f\n',sqrt(rhon1(k)));
   
end
if k>1 
    beta=sqrt(rhon1(k-1));
else
    beta=sqrt(rhon1(1));
end
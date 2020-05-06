function [p]=proj_level(x,eta,beta)
if eta>=1
    p=x;
elseif norm(x)<=sqrt(eta/beta)
    p=x;
else
    p=x/norm(x)*sqrt(eta/beta);
end
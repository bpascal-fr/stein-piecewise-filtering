function wp = prox_L1(wx,gamma)
% Proximity operator of the l1 norm 
% wp = prox_{gamma || .||_1}(wx)
%
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% June 2019

wp = max(abs(wx)-gamma,0).*sign(wx);

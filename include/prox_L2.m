function wp = prox_L2(wx,gamma)
% Proximity operator of the l2 norm 
% wp = prox_{gamma || .||_2^2}(wx)
%
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% June 2019

wp = wx/(1+gamma);

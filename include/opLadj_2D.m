function x = opLadj_2D(y1,y2)
% Adjoint of the operator associated with horizontal and vertical
% difference operator.
%
% Implementation N. PUSTELNIK, CNRS, ENS Lyon
% June 2019

tau = 1;
y = y1;
[n,m] = size(y);
x = -([y(:,1:tau)/2,y(:,1+tau:m-tau)/2- y(:,1:m-2*tau)/2,-y(:,m-2*tau+1:m-tau)/2]);
y = y2';
[n,m] = size(y);
x = x -([y(:,1:tau)/2,y(:,1+tau:m-tau)/2- y(:,1:m-2*tau)/2,-y(:,m-2*tau+1:m-tau)/2])';


function [F] = root4d(x,lambda,k,gamma,phi, nonlinearity)
% 
%
% F(1) = Re{A} = x(1)
% F(2) = Re{B} = x(2)
% F(3) = Im{A} = x(3)
% F(4) = Im{B} = x(4)
%
% k - coupling strength    
% gammas - gain and loss vector of parameters gamma.
% phi - Peierls phase
% nonlinearity - nonlinearity strength
%
% Author: Claudia Castro-Castro @ SMU
% Date  : September 2016

% solving for real parts
F(1) = -gamma*x(3) + (-lambda + nonlinearity*(x(1)^2 + x(3)^2) )*x(1) + 2*k*cos(phi)*x(2);
F(2) =  gamma*x(4) + (-lambda + nonlinearity*(x(2)^2 + x(4)^2) )*x(2) + 2*k*cos(phi)*x(1);


% solving for imaginary parts
F(3) =  gamma*x(1) + (-lambda + nonlinearity*(x(1)^2 + x(3)^2) )*x(3) + 2*k*cos(phi)*x(4);
F(4) = -gamma*x(2) + (-lambda + nonlinearity*(x(2)^2 + x(4)^2) )*x(4) + 2*k*cos(phi)*x(3);


end
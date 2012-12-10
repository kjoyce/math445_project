function xhat = ABCexact(time,theta,y0)
%%%%%%%%%%%%%%
% This is not the correct model.  You need to solve for Xa in the
% differential equation
%%%%%%%%%%%%%
    y0 = y0(1);
    xhat = theta(1)*y0/(theta(1) - theta(2))*(exp(-theta(2)*time) - exp(-theta(1)*time));
    xhat = xhat';
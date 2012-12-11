function xhat = ABCexact(time,theta,y0)
%%%%%%%%%%%%%%
% This is not the correct model.  You need to solve for Xa in the
% differential equation
%%%%%%%%%%%%%
    y0 = y0(1);
    time = time'; % We need to use columns now!
    xhat = zeros(length(time),2);
    %%% CHANGE TO ACTUAL SOLUTION WHEN YOU FIND IT
    xhat(:,1) = theta(3)*y0*exp(-theta(2)*time);
    xhat(:,2) = theta(1)*y0/(theta(1) - theta(2))*(exp(-theta(2)*time) - exp(-theta(1)*time));
 

function dy=ABCode(t,y,theta)

% take parameters and components out from y and theta
k1=theta(1); k2=theta(2);
A=y(1); B=y(2);

% define the ODE
dy(1) = -k1*A;
dy(2) = k1*A-k2*B;
dy(3) = k2*B;

dy=dy(:);  % make sure that we return a column vector
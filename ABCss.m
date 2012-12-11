function ss=ABCss(theta,mouse,varargin)

time=[0 mouse.xdata];  % adding zero to time vector!
ydata=mouse.ydata;
y0=mouse.init;

if (length(varargin{:})>1) ...
&& strcmp('use_ode_solution',varargin{1}{1})...
&& (varargin{1}{2} == 1)
    ymod=ABCexact(time,theta,y0);    
else
    ymod=ABCmodel(time,theta,y0);
end

ymod=theta(3)*ymod(2:end,2);     % Fit only the X (blood) data points

ss=sum(sum((ydata-ymod').^2));   % the total SS
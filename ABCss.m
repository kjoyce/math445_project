function ss=ABCss(theta,data)

time=[0 data.time];  % adding zero to time vector!
ydata=data.ydata;
y0=data.y0;

ymod=ABCmodel(time,theta,y0);
ymod=theta(3)*ymod(2:end,2);  % taking A and B, removing initial value

ss=sum(sum((ydata-ymod).^2));   % the total SS
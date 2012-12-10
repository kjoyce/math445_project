function ss = monodss(theta,xdata,ydata);

y = theta(1)*xdata./(theta(2)+xdata);
ss = sum((y-ydata).^2);

function ss = ABCsumss(theta,mice)
ss=0;
 t = linspace(0,10);
 plot(t,ABCmodel(t,theta,[mice(1).init 0]),mice(1).xdata,mice(1).ydata,'ko');
 drawnow;
for i=1:length(mice)
    data.time = mice(i).xdata;
    data.ydata = mice(i).ydata';
    data.y0 = [mice(i).init 0];
    
    % Note we are not using oldpar
    ss     = ss + ABCss(theta,data); % first SS value
end

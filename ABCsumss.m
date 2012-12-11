function ss = ABCsumss(theta,mice,varargin)
ss=0;
%  t = linspace(0,10);
%  plot(t,ABCmodel(t,theta,[mice(1).init 0]),mice(1).xdata,mice(1).ydata,'ko');
%  drawnow;
for i=1:length(mice)
    ss     = ss + ABCss(theta,mice(i),varargin); % first SS value
end

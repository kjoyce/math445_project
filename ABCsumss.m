function ss = ABCsumss(theta,mice)
ss=0;
for i=1:length(mice)
    data.time = mice(i).xdata;
    data.ydata = mice(i).ydata';
    data.y0 = [mice(i).init 0];
    
    % Note we are not using oldpar
    ss     = ss + ABCss(theta,data); % first SS value
end

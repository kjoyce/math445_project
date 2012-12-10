%clear all, close all, clc;

load mice_data_struct;
the_mice = mice.b;
figure();
hold on;
for (i=1:10)
    data.time = the_mice(i).xdata;
    data.ydata = the_mice(i).ydata';
    data.y0 = [the_mice(i).init 0];
    
    % calling fminsearch via some magic from MCMC
    theta0=[1 1];
    [thopts(i,:),ssopt]=fminsearch(@ABCexactss,theta0,[],data);
    
    % visualization: solve model with thopt and compare to data
    t=linspace(0,10);
    ymod=ABCexact(t,thopts(i,:),data.y0);
    
    % figure()
    plot(t,ymod);% hold on;
    plot(data.time,data.ydata,'o'); %hold off;
    xlabel('time'); ylabel('Drug Concentration');
    %legend('Xa compartment','X Compartment');
    % hold off;
end
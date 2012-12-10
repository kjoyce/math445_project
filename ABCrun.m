%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK: Fit cascading absorbtion model to mice data.
% 
% Model:
% 
% dXa/dt = 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A --> B --> C demo
%clear all; close all; clc;

load mice_data_struct;
% the data structure
%%% CHANGE INTO LOADING OUR DATA
% for (i,1:20)
% data.time  = 1:2:9;
% data.ydata = [.504	.415
%               .217	.594
%               .101	.493
%               .064	.394
%               .008   .309];
% data.y0    = [1 0 0];

data.time = mice.a(1).xdata;
data.ydata = mice.a(1).ydata';
data.y0 = [mice.a(1).init 0];

% calling fminsearch
theta0=[1 1];
[thopt,ssopt]=fminsearch(@ABCss,theta0,[],data);

% visualization: solve model with thopt and compare to data
t=linspace(0,10);
ymod=ABCmodel(t,thopt,data.y0);

plot(t,ymod); hold on;
plot(data.time,data.ydata,'o'); hold off;
xlabel('time'); ylabel('Drug Concentration');
legend('Xa compartment','X Compartment');

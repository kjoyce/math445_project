function y=ABCmodel(time,theta,y0)

[t,y]=ode45(@ABCode,time,y0,[],theta);
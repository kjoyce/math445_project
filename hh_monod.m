

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK: Repeat this run with the following example,  from 
% P. M. Berthouex and L. C. Brown:
% _Statistics for Environmental Engineers_, CRC Press, 2002.
%
% Fit the Monod model
%
%  y = theta_1 * x /(theta_2 + x)  + eps,  epsilon ~ N(0,I\sigma^2) 
%
% to observations
%
%   x (mg / L COD):  28    55    83    110   138   225   375   
%   y (1 / h):       0.053 0.060 0.112 0.105 0.099 0.122 0.125
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% The task:
xdata = [ 28    55    83    110   138   225   375  ];
ydata = [ 0.053 0.060 0.112 0.105 0.099 0.122 0.125];


%first plot the data,  
plot(xdata,ydata,'o-');
%break

% the value of the Monod model for large 'x' is close to 'theta_1', so from
% the plot get the initial guess theta_1=0.1 , for theta_2 just take something 
% of the size of 'xdata'. So take the initial guess by hte plot:

par_ini = [.1 .1];    % givew here the initial guess for fit

nsimu  = 5000;                     % length of chain
% call of the optimizer:
[par0,rss]   = fminsearch('monodss',par_ini,[],xdata,ydata);
plot(xdata,ydata,'o', xdata,par0(1)*xdata./(par0(2)+xdata),'r-');
%pause; 

% get the estimate for sigma2:
nobs = length(xdata);
npar = 2;                            %  n of parameters,th1 ja th2
sigma2 = rss/(nobs-npar);            %  variance of meas.error

%As for choice of the proposal, try the following:

%the 'initial guess' for the proposal:
qcov   = 1e-4*eye(npar);              %  covariance for Gaussian proposal
%qcov = cov(chain);

% ... but when one (or more) run has been done, comment the above and try
% to find a better covariance for a Gaussian proposal

%%%%%%%%%%%%%%%MCMC sampling starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R      = chol(qcov);
oldpar = par0;
ss     = monodss(oldpar,xdata,ydata); % first SS value
rej        = 0;                     % initialized count for rejections

chain(1,:) = oldpar;
for i=2:nsimu % simulation loop
   newpar = oldpar+randn(1,npar)*R;     % new candidate
   ss2    = ss;                         % old SS
   ss1    = monodss(newpar,xdata,ydata);  % new SS
   if (ss1<ss2 | rand < exp(-0.5*(ss1-ss2)/sigma2))
     chain(i,:) = newpar;              % accept
     oldpar     = newpar;
     ss         = ss1;   
   else   
     chain(i,:) = oldpar;              % reject
     rej        = rej+1;  
   end
end

accept = 1-rej./nsimu;                  % acceptance rate

%%%%%%%%%%%MCMC sampling done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%some plots: 
figure(1)
subplot(2,1,1);plot(1:nsimu,chain(:,1)); title('CHAIN for theta_1');
subplot(2,1,2);plot(1:nsimu,chain(:,2)); title('CHAIN for theta_2');

figure(2); plot(chain(:,1),chain(:,2),'.');
title('SAMPLED PARAMETER POSTERIOR DISTRIBUTION');


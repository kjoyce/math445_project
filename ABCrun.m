%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK: Fit cascading absorbtion model to mice data.
% 
% Model:
% 
% dXa/dt = 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A --> B --> C demo
clear all; close all; clc;

load mice_data_struct;

% CHANGE THESE PARAMETERS %
nsimu = 5000
the_mice = mice.a;     

% calling fminsearch
theta0=[10 1 1];
[thopt,ssopt]=fminsearch(@ABCsumss,theta0,[],the_mice);

% visualization: solve model with thopt and compare to data
t=linspace(0,10);
figure();
hold on;
for i=1:10
    ymod=ABCmodel(t,thopt,[the_mice(i).init 0]);
    
    plot(t,ymod);% hold on;
    plot(the_mice(i).xdata,the_mice(i).ydata,'ko'); %hold off;
    xlabel('time'); ylabel('Drug Concentration');
    %legend('Xa compartment','X Compartment');
    % hold off;
end

ss = ssopt;
oldpar = thopt;
nobs = 40;
npar = 3;                         %  n of parameters,th1 ja th2
sigma2 = ss/(nobs-npar);          %  variance of meas.error

%As for choice of the proposal, try the following:

%the 'initial guess' for the proposal:
%qcov   = 1e-4*eye(npar);              %  covariance for Gaussian proposal
%qcov = cov(chain);

%... but when one (or more) run has been done, comment the above and try
%to find a better covariance for a Gaussian proposal

%%%%%%%%%%%%%%MCMC sampling starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R      = chol(qcov);
rej        = 0;                     % initialized count for rejections

chain = zeros(nsimu,npar);
chain(1,:) = oldpar;
for i=2:nsimu % simulation loop  
  newpar = oldpar+randn(1,npar)*R;     % new candidate
  ss2    = ss;                         % old SS
  ss1    = ABCsumss(newpar,the_mice);  % new SS
  if (ss1<ss2 | rand < exp(-0.5*(ss1-ss2)/sigma2))
    chain(i,:) = newpar;              % accept
    oldpar     = newpar;
    ss         = ss1;   
  else   
    chain(i,:) = oldpar;              % reject
    rej        = rej+1;  
  end
  if (mod(i,10) == 0)
      fmt = 'MCMC run %i\n';
      fprintf(fmt,i);
  end
  %if (i > 100) 
  if( mod(i,50) == 0 ) %speeds up the MCMC
    qcov = cov(chain);
    R = chol(qcov);
  end
end

accept = 1-rej./nsimu;                  % acceptance rate

%%%%%%%%%%%MCMC sampling done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%some plots: 
figure()
subplot(3,1,1);plot(1:nsimu,chain(:,1)); title('CHAIN for ka');
subplot(3,1,2);plot(1:nsimu,chain(:,2)); title('CHAIN for K');
subplot(3,1,3);plot(1:nsimu,chain(:,3)); title('CHAIN for F/V');


figure(); 
subplot(2,2,1);plot(chain(:,1),chain(:,2),'.'); title('ka vs K');
subplot(2,2,2);plot(chain(:,1),chain(:,3),'.'); title('ka vs F/V');
subplot(2,2,3);plot(chain(:,2),chain(:,3),'.'); title('K vs F/V');
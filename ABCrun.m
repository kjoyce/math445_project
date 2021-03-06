%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK: Fit cascading absorbtion model to mice data.
% 
% Model:
% 
% dXa/dt = -ka * Xa
% dX/dt  = ka * Xa - K * X
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leave this commented, unless you need to do some debugging
%clear all; close all; clc; 

%%%%%%%%%%%% CHANGE THESE PARAMETERS %%%%%%%%%%%%%%%%
load mice_data_struct;
nsimu = 15000;
mice_group = 'a';
%mice_group = 'b';
first_run = 1;  % change to 0 if you are continuing from previous MCMC
%%% If first_run = 0, then define the chain_file_name 
% mice_a_15000-20000_exact-sol.mat
% You probably won't ever have to do this since the chain calculation is blazing fast now.
if strcmp(mice_group,'a')
  the_mice = mice.a;
  chain_file_name = ['mice_a_1-',num2str(nsimu),'_exact-sol.mat'];
  fig_save_prefix = 'mice_a_';
  plot_label = 'Mice A'
elseif strcmp(mice_group,'b')
  the_mice = mice.b;
  chain_file_name = ['mice_b_1-',num2str(nsimu),'_exact-sol.mat'];
  fig_save_prefix = 'mice_b_';
  plot_label = 'Mice B'
else
  error('Bad Mice Choice');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if (first_run)
  % calling fminsearch
  theta0=[10 1 1];
  [thopt,ssopt]=fminsearch(@ABCsumss,theta0,[],the_mice);

  % visualization: solve model with thopt and compare to data
  t=linspace(0,10);
  figure(1);
  set(gcf,'PaperSize',[5,4]);
  set(gcf,'PaperPosition',[0,0,5,4]);
  hold on;
  for i=1:10
      ymod=ABCmodel(t,thopt,the_mice(i).init);
    %  plot(t,ymod(:,1),':c');
    %  plot(t,ymod(:,2),':r');
      plot(t,ymod);
      plot(the_mice(i).xdata,the_mice(i).ydata,'ko');
      xlabel('time'); ylabel('Drug Concentration');
  end
%      legend('Xa compartment','X Compartment');
  hold off;

  figure(5);
  set(gcf,'PaperSize',[5,4]);
  set(gcf,'PaperPosition',[0,0,5,4]);
  idx = randperm(10); % plot nine random mice
  for i=1:9
      ymod=ABCmodel(t,thopt,the_mice(idx(i)).init);
      subplot(3,3,i);
      hold on;
      plot(t,ymod);
      plot(the_mice(i).xdata,the_mice(idx(i)).ydata,'ko');
      ylim([-1,40]);
      hold off;
      xlabel('time'); ylabel('Drug Concentration');
  end
  print(gcf,'-dpdf',[fig_save_prefix 'fits.pdf']) 
end    
nobs = 40;
npar = 3;                         %  n of parameters,th1 ja th2
if (first_run)
    %As for choice of the proposal, try the following:
    %the 'initial guess' for the proposal:
    qcov   = 1e-4*eye(npar);          %  covariance for Gaussian proposal
    ss = ssopt;
    sigma2 = ss/(nobs-npar);          %  variance of meas.error
    oldpar = thopt;
    chain = zeros(nsimu,npar);        % initialize the memory for chain (slight speedup)
    chain(1,:) = oldpar;              % initialize the chain
    start_chain = 1;
else % When this is run be sure that you've loaded ss, chain, oldpar, and qcov
    start_chain = length(chain)+1;
    chain = [chain; zeros(nsimu,npar)];% initialize the memory for chain (slight speedup)
end                                     
                                      
%%%%%%%%%%%%%%MCMC sampling starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R      = chol(qcov);
rej    = 0;                         % initialized count for rejections
for i=start_chain:start_chain+nsimu % simulation loop  
  newpar = oldpar+randn(1,npar)*R;     % new candidate
  ss2    = ss;                         % old SS
  ss1    = ABCsumss(newpar,the_mice,'use_ode_solution',1);  % new SS
  if (ss1<ss2 || rand < exp(-0.5*(ss1-ss2)/sigma2))
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
  if (i > 100) 
  %if( mod(i,50) == 0 ) % Less adapting, speeds up when chain is superlong
    qcov = cov(chain(1:i,:));
    R = chol(qcov);
  end
end

accept = 1-rej./nsimu;                  % acceptance rate

%%%%%%%%%%%MCMC sampling done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%some plots: 
figure(3)
set(gcf,'PaperSize',[5,4]);
set(gcf,'PaperPosition',[0,0,5,4]);
start_chain = 1;
chlen = length(chain);
subplot(3,1,1);plot(start_chain:chlen,chain(start_chain:end,1)); title('CHAIN for ka'); xlim([0,chlen]);
ylim([min(chain(:,1)),max(chain(:,1))]);
subplot(3,1,2);plot(start_chain:chlen,chain(start_chain:end,2)); title('CHAIN for K'); xlim([0,chlen]);  
ylim([min(chain(:,2)),max(chain(:,2))]);
subplot(3,1,3);plot(start_chain:chlen,chain(start_chain:end,3)); title('CHAIN for F/V'); xlim([0,chlen]);%  ylim([.9,1.2]);
ylim([min(chain(:,3)),max(chain(:,3))]);
print(gcf,'-dpdf',[fig_save_prefix 'chains.pdf']) 

figure(4); 
set(gcf,'PaperSize',[5,4]);
set(gcf,'PaperPosition',[0,0,5,4]);
subplot(2,2,1);plot(chain(start_chain:end,1),chain(start_chain:end,2),'.'); title('ka vs K'); 
subplot(2,2,2);plot(chain(start_chain:end,1),chain(start_chain:end,3),'.'); title('ka vs F/V');
subplot(2,2,3);plot(chain(start_chain:end,2),chain(start_chain:end,3),'.'); title('K vs F/V');
subplot(2,2,4); axis off;
text(.1,.5,plot_label,'FontSize',32);
print(gcf,'-dpdf',[fig_save_prefix 'joint_dists.pdf']) 


save(chain_file_name);


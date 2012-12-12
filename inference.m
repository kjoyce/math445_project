load 'mice_a_1-15000_exact-sol.mat';
chain_a = chain;
load 'mice_b_1-15000_exact-sol.mat';
chain_b = chain;

diff_chain = chain_a - chain_b;

diff_K = diff_chain(:,2);
diff_FV = diff_chain(:,3);

figure(1)
set(gcf,'PaperSize',[5,4]);
set(gcf,'PaperPosition',[0,0,5,4]);
hist(diff_K,100)
title('Parameter Distribution of The Drug Elimination Rate Constant Difference')

figure(2)
set(gcf,'PaperSize',[5,4]);
set(gcf,'PaperPosition',[0,0,5,4]);
hist(diff_FV,100)
title('Parameter Distribution of The Drug Absorption Difference')

alpha_lev = .05;
l_idx = floor(alpha_lev/2*length(diff_chain));
r_idx = ceil((1-alpha_lev/2)*length(diff_chain));
sort_diff_K = sort(diff_K);
CI95_diff_K = [sort_diff_K(l_idx) sort_diff_K(r_idx)] 
sort_diff_FV = sort(diff_FV);
CI95_diff_FV = [sort_diff_FV(l_idx) sort_diff_FV(r_idx)] 

% Checking for autocorrelation
figure(3);
set(gcf,'PaperSize',[5,4]);
set(gcf,'PaperPosition',[0,0,5,4]);
hist(diff_K(1:10:end),50);
title('Every 10th entry of the difference chain for K');
print(gcf,'-dpdf','auto_cor_K.pdf')

figure(4);
set(gcf,'PaperSize',[5,4]);
set(gcf,'PaperPosition',[0,0,5,4]);
hist(diff_FV(1:10:end),50);
title('Every 10th entry of the difference chain for F/V');
print(gcf,'-dpdf','auto_cor_FV.pdf')

alpha_lev = .05;
l_idx = floor(alpha_lev/2*length(diff_chain)/10);
r_idx = ceil((1-alpha_lev/2)*length(diff_chain)/10);
sort_diff_K = sort(diff_K(1:10:end));
CI95_diff_K = [sort_diff_K(l_idx) sort_diff_K(r_idx)] 
sort_diff_FV = sort(diff_FV(1:10:end));
CI95_diff_FV = [sort_diff_FV(l_idx) sort_diff_FV(r_idx)] 


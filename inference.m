load 'mice_a_1-15000_exact-sol.mat';
chain_a = chain;
load 'mice_b_1-15000_exact-sol.mat';
chain_b = chain;

diff_chain = chain_a - chain_b;

diff_K = diff_chain(:,2);
diff_FV = diff_chain(:,3);

hist(diff_K,50)
title('Parameter Distribution of The Drug Elimination Rate Constant Difference')

hist(diff_FV,50)
title('Parameter Distribution of The Drug Absorption Difference')

alpha_lev = .05;
l_idx = ceil(alpha_lev/2*length(diff_chain));
r_idx = ceil((1-alpha_lev/2)*length(diff_chain));
sort_diff_K = sort(diff_K);
CI95_diff_K = [sort_diff_K(l_idx) sort_diff_K(r_idx)] 
sort_diff_FV = sort(diff_FV);
CI95_diff_FV = [sort_diff_FV(l_idx) sort_diff_FV(r_idx)] 

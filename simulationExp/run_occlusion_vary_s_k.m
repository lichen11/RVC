% Assessing Robustness
% occlusion contamination
addpath '~/Dropbox/Research/General_code/';
p = 0.9;
save_file_name = sprintf('occlusion_pt%d_vary_s_k_d4.mat', 10*p);
 
rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
n = 200;
%d_vec = 1:20;
k_vec = 1:20;
s_vec = k_vec;
nsim= 500;
occlusionErrTensor = zeros(length(s_vec), 3, nsim);

parfor i = 1:nsim
    rng(i);
    occlusionErr = zeros(length(s_vec), 3);
    [A, tau, ~] = sbm(rho, B, n);
    for j = 1: length(k_vec)
        k = k_vec(j);
        s = k;        
        occlusionErr_ = occlusion(p, A, tau, 4, s, k);
        occlusionErr(j,:) = occlusionErr_;
    end
    occlusionErrTensor(:,:,i) = occlusionErr;
end

save(save_file_name, 'rho', 'B', 'n', 'p', 'nsim', 'k_vec', 's_vec', 'occlusionErrTensor');
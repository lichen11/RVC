% Assessing Robustness
% occlusion contamination
%addpath '~/Dropbox/Research/General_code/';
addpath '../Real_data_exp/';
p = 0.8;
save_file_name = sprintf('run_occlusion_adjnoun_like_pt%d_vary_s_d.mat', 10*p);
load estimateB_adjnoun
%rho = [0.4 0.6];
%B = [0.7 0.32; 0.32 0.75]; 
n = 200;
d_vec = 1:60;
%k_vec = 1:20;
K =1;
s_vec = d_vec;
nsim= 100;
occlusionErrTensor = zeros(length(s_vec), 3, nsim);

parfor i = 1:nsim
    rng(i);
    occlusionErr = zeros(length(s_vec), 3);
    [A, tau, ~] = sbm(rho, B, n);
    for j = 1: length(d_vec)
        d = d_vec(j);
        s = d;        
        occlusionErr_ = occlusion(p, A, tau, d, s, K);
        occlusionErr(j,:) = occlusionErr_;
    end
    occlusionErrTensor(:,:,i) = occlusionErr;
end

save(save_file_name, 'rho', 'B', 'n', 'p', 'nsim', 'd_vec', 's_vec', 'occlusionErrTensor');
% Assessing Robustness
% occlusion reversion contamination
addpath '~/Dropbox/Research/General_code/';
p = 0.8;
save_file_name = sprintf('run_occlusion_reversion_pt%d_vary_s_d.mat', 10*p);

rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
n = 200;
d_vec = 1:20;
s_vec = d_vec;
nsim= 500;
occlusion_reversion_ErrTensor = zeros(length(d_vec), 3, nsim);

parfor i = 1:nsim
    rng(i);
    occlusion_reversion_Err = zeros(length(d_vec), 3);
    [A, tau, ~] = sbm(rho, B, n);
    for j = 1: length(d_vec)
        d = d_vec(j);
        s = d;        
        occlusion_reversion_Err_ = mixContamination(p, A, tau, d, s, 'occlusion_reversion');
        occlusion_reversion_Err(j,:) = occlusion_reversion_Err_;
    end
    occlusion_reversion_ErrTensor(:,:,i) = occlusion_reversion_Err;
end

save(save_file_name, 'rho', 'B', 'n', 'p', 'nsim', 'd_vec', 's_vec', 'occlusion_reversion_ErrTensor');
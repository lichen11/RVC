% Assessing Robustness
% occlusion reversion contamination
addpath '~/Dropbox/Research/General_code/';
%p = 0.8;
save_file_name = 'run_occlusion_reversion.mat';%sprintf('run_occlusion_reversion_pt%d_vary_s_d.mat', 10*p);

rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
n = 200;
p_vec = 0: 0.03: 0.93;

nsim= 500;
occlusion_reversion_ErrTensor1 = zeros(length(p_vec), 3, nsim);
occlusion_reversion_ErrTensor2 = zeros(length(p_vec), 3, nsim);

parfor i = 1:nsim
    rng(i);
    occlusion_reversion_Err1 = zeros(length(p_vec), 3);
    occlusion_reversion_Err2 = zeros(length(p_vec), 3);
    [A, tau, ~] = sbm(rho, B, n);
    for j = 1: length(p_vec)
        p = p_vec(j);        
        occlusion_reversion_Err1_ = mixContamination(p, A, tau, 2, 5, 'occlusion_reversion');
        occlusion_reversion_Err1(j,:) = occlusion_reversion_Err1_;
        occlusion_reversion_Err2_ = mixContamination(p, A, tau, 4, 5, 'occlusion_reversion');
        occlusion_reversion_Err2(j,:) = occlusion_reversion_Err2_;
    end
    occlusion_reversion_ErrTensor1(:,:,i) = occlusion_reversion_Err1;
    occlusion_reversion_ErrTensor2(:,:,i) = occlusion_reversion_Err2;

end

save(save_file_name, 'rho', 'B', 'n', 'nsim', 'p_vec', 'occlusion_reversion_ErrTensor1', 'occlusion_reversion_ErrTensor2');
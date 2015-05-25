addpath '~/Dropbox/Research/General_code/';
p = 0.6;
save_file_name = sprintf('reversion_pt%d_vary_s_d_newB.mat', 10*p);
rho = [0.4 0.6];
%B = [0.7 0.32; 0.32 0.75]; 
B = [0.2 0.1; 0.1 0.35]; 
n = 200;
nsim = 500;
d_vec = 1:20;
s_vec = d_vec;
flipErrTensor = zeros(length(d_vec), 3, nsim);

parfor nexp = 1 : nsim
    rng(nexp);
    flipErr = zeros(length(d_vec), 3);
    [A, tau, ~] = sbm(rho, B, n);
    for j = 1: length(s_vec)
        d = d_vec(j);
        s = s_vec(j);
        flipErr_ = flipping(p, A, tau, d, s);
        flipErr(j,:) = flipErr_;        
    end
    flipErrTensor(:,:,nexp) = flipErr;
end

save(save_file_name, 'rho', 'B', 'n', 'p', 'nsim', 'd_vec', ...
    's_vec', 'flipErrTensor');
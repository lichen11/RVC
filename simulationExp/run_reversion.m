%rho = [0.4 0.6];
%B = [0.7 0.32; 0.32 0.75]; 
load estimateB_adjnoun
n = 200;
p = 0: 0.03: 1;
nsim = 100;
save_file_name = 'run_reversion_sparsity_embedding_adjnoun_like.mat';
flipErrTensor = zeros(length(p), 3, nsim);
flipErrTensor4D = zeros(length(p), 3, nsim);
parfor nexp = 1 : nsim
    rng(nexp);
    [A, tau, ~] = sbm(rho, B, n);
    flipErr = flipping(p, A, tau, 2, 5);
    flipErrTensor(:, :, nexp) = flipErr;
    flipErr4D = flipping(p, A, tau, 4, 6);
    flipErrTensor4D(:, :, nexp) = flipErr4D;
end

save(save_file_name, 'B', 'rho', 'n', 'p','nsim','flipErrTensor','flipErrTensor4D');

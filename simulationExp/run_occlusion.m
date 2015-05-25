% Assessing Robustness
% occlusion contamination
%rho = [0.4 0.6];
%B = [0.7 0.32; 0.32 0.75];
load estimateB_adjnoun
n = 200;
p = 0: 0.03: 0.93;
save_file_name = 'run_occlusion_sparsity_embedding_adjnoun_like.mat';
nsim= 100;
occlusionErrTensor = zeros(length(p), 3, nsim);
occlusionErrTensorSparsity3 = zeros(length(p), 3, nsim);
occlusionErrTensorSparsity4 = zeros(length(p), 3, nsim);
occlusionErrTensorSparsity5 = zeros(length(p), 3, nsim);
occlusionErrTensor4D = zeros(length(p), 3, nsim);
parfor i = 1:nsim
    rng(i);
    [A, tau, ~] = sbm(rho, B, n);
    occlusionErr = occlusion(p, A, tau, 2, 5, 1);
    %occlusionErrSparsity3 = occlusion(p, A, tau, 2, 3, 1);
    %occlusionErrSparsity4 = occlusion(p, A, tau, 2, 4, 1);
    %occlusionErrSparsity5 = occlusion(p, A, tau, 2, 5, 1);
    occlusionErr4D = occlusion(p, A, tau, 4, 5, 1);
    occlusionErrTensor(:,:,i) = occlusionErr;
    %occlusionErrTensorSparsity3(:,:,i) = occlusionErrSparsity3;
    %occlusionErrTensorSparsity4(:,:,i) = occlusionErrSparsity4;
    %occlusionErrTensorSparsity5(:,:,i) = occlusionErrSparsity5;
    occlusionErrTensor4D(:,:,i) = occlusionErr4D;
end

%save(save_file_name, 'rho', 'B', 'n', 'p', 'nsim', 'occlusionErrTensor', 'occlusionErrTensorSparsity3',...
%    'occlusionErrTensorSparsity4','occlusionErrTensorSparsity5', 'occlusionErrTensor4D');
save(save_file_name, 'rho', 'B', 'n', 'p', 'nsim', 'occlusionErrTensor', 'occlusionErrTensor4D');
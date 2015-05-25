addpath '~/Dropbox/Research/General_code';
addpath '~/Dropbox/Research/Robust_est/RealDataExpNew/jitter';
addpath '~/Dropbox/Research/Robust_est/Real_data_exp/enron';
addpath(genpath('~/Dropbox/FamousNetworks/.'))
addpath(genpath('../simulationExp/.'))
% load data file
load '~/Dropbox/VN/data/celegans_no_isolates_symmetric_binarized.mat';
%load 'Label'
%save file name
save_file_name = 'run_SRC_on_celegans_1_100.mat';
A = A;%Adj;
Lab = labels;
% sparsity variation
s_vec = 1:100;
d_vec = 1:100;
errMat = zeros(length(d_vec), 5); % SRC, 1NN, 3NN, 5NN, LDA

for j = 1:length(d_vec)
    d = d_vec(j);
    s = s_vec(j);
    % SRC for various sparsity
    tmp = errfun(Lab, A, s);  
    errMat(j,1) = tmp;
    % diagonal augmentation and ASE, apply 1 NN
    ase = adjSpecEmbed(diagAug( A ), d);
    tmpknn1 = errsvdknn(Lab, ase, 1); 
    tmpknn3 = errsvdknn(Lab, ase, 3); 
    tmpknn5 = errsvdknn(Lab, ase, 5); 
    errMat(j, 2:4) =  [tmpknn1 tmpknn3 tmpknn5];
    % apply LDA
    det_tmp = det(cov(ase));
    if det_tmp < 10^(-10)
        tmplda =  errsvdlda(Lab, awgn(ase, 250));
    else
        tmplda =  errsvdlda(Lab, ase);
    end
    errMat(j, 5) = tmplda;
end

save(save_file_name, 'errMat', 's_vec', 'd_vec')


%figure; plot( 1:dlim, srcErr, '-*b', 1:dlim,  smooth(errVaryingD(1:dlim, 1)), '-ok',...
%    1:dlim, smooth(errVaryingD(1:dlim,2)),  '-^r', 'LineWidth', 2, 'MarkerSize',7)

%legend('SRC','NN o ASE','LDA o ASE','Location','SouthWest');
%set(gca,'FontSize',18)
%xlabel('Embedding Dimension','fontsize',32)
%ylabel('Error rate','fontsize',32)
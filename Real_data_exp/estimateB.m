% estimate the B matrix from real data
addpath '~/Dropbox/Research/General_code';
addpath '~/Dropbox/Research/Robust_est/RealDataExpNew/jitter';
addpath '~/Dropbox/Research/Robust_est/Real_data_exp/enron';
addpath(genpath('~/Dropbox/FamousNetworks/.'))

load 'adjnoun.mat';
%load 'Label'
%save file name
save_file_name = '~/Dropbox/Research/Robust_est/Real_data_exp/estimateB_adjnoun.mat';
A = Adj;%Adj;
Lab = Label+1;


nClass = length(unique(Lab));
tb = tabulate(Lab);
m = tb(:,2);
rho = tb(:,3)/100;
B= zeros(nClass);
for k = 1:nClass
    for j = 1: nClass
        B(k, j) = sum(sum(A(Lab==k, Lab==j)))/m(k)/m(j);
    end  
end

for k = 1:nClass
        B(k, k) = m(k)/(m(k) - 1)*B(k, k);
end

save(save_file_name, 'B', 'rho')

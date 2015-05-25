% get the table that says this is a good paper!
addpath '~/Dropbox/Research/General_code';
%nBlock = 2; 
warning('off','all')
numsim = 300;
rho = [0.6 0.4];
B = [0.42 0.42; 0.42 0.65]; 
n_vec = 20:10:250;
err = zeros(length(n_vec), 3, numsim);
s = 5;
K = 1;
d = 2;
save_file_name = 'simulation_no_contamination_same_cov2.mat';
parfor i = 1:numsim
	err_ = zeros(length(n_vec), 3);
	rng(i)
	for j = 1:length(n_vec)
		n = n_vec(j);
		[A, tau, ~] = sbm(rho, B, n);
		%A = diagAug(A);
		embed = svdembed(d, A);
		exp = 1;
		% SRC
		[err_(j, exp), ~] = errfun(tau, A, s);
		exp = exp + 1;		
		% KNN		
		[err_(j, exp), ~] = errsvdknn(tau, embed, K);
		exp = exp + 1;	
		% LDA		
		[err_(j, exp), ~] = errsvdlda(tau, embed);
	end
	err(:,:,i) = err_;
end
save(save_file_name, 'err', 'rho', 'B', 'n_vec', 'numsim', 'K', 's');

% run synthetic data based on real data performance

addpath '~/Dropbox/Research/General_code';
addpath '../Real_data_exp/';
load estimateB_adjnoun
%nBlock = 2; 
warning('off','all')
numsim = 500;
%rho = [0.4 0.6];
%B = [0.0091 0.0002; 0.2 0.75]; 
n = 120;
d_vec = 1:3:90;
K = 1;
s_vec = d_vec;
err = zeros(length(d_vec), 3, numsim);
save_file_name = 'simulation_adjnounlike_n_120_no_contamination_s_d.mat';
parfor i = 1:numsim
	err_ = zeros(length(d_vec), 3);
	rng(i)
	for j = 1:length(d_vec)
		d = d_vec(j);		
		s = d_vec(j);
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
        det_tmp = det(cov(embed(:,1:d)));
        if det_tmp < 10^(-10)
            tmplda =  errsvdlda(tau, awgn( embed(:,1:d), 250));
        else
            tmplda =  errsvdlda(tau,  embed(:,1:d));
        end
		err_(j, exp) = tmplda;
		%[err_(j, exp), ~] = errsvdlda(tau, embed);
	end
	err(:,:,i) = err_;
end
save(save_file_name, 'err', 'rho', 'B', 'n', 'numsim', 'K', 's_vec', 'd_vec');

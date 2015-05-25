% get the table that says this is a good paper!
addpath '~/Dropbox/Research/General_code';
%nBlock = 2; 
warning('off','all')
numsim = 200;%500;
rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
n = 110;
d_vec = 1:3:52;
K = 1;
s_vec = 1:3:52;
err_ase = zeros(length(d_vec), 2, numsim);
err_src = zeros(length(s_vec), numsim);
save_file_name = 'simulation_n_100_no_contamination_s_d.mat';
parfor i = 1:numsim
	err_ase_ = zeros(length(d_vec), 2);
    err_src_ = zeros(length(s_vec), 1);
	%rng(i)
    [A, tau, ~] = sbm(rho, B, n);
    d = d_vec(end);
    embed = svdembed(d, A);
    
    for j = 1:length(s_vec)
        s = s_vec(j);
		[err_src_(j), ~] = errfun(tau, A, s);
    end
	for j = 1:length(d_vec)
		d = d_vec(j);		
		%A = diagAug(A);
		exp = 1;	
		% KNN		
		[err_ase_(j, exp), ~] = errsvdknn(tau, embed(:,1:d), K);
		exp = exp + 1;	
		% LDA		
        det_tmp = det(cov(embed(:,1:d)));
        if det_tmp < 10^(-10)
            tmplda =  errsvdlda(tau, awgn( embed(:,1:d), 250));
        else
            tmplda =  errsvdlda(tau,  embed(:,1:d));
        end
		err_ase_(j, exp) = tmplda; %errsvdlda(tau, embed(:,1:d));
	end
	err_ase(:,:,i) = err_ase_;
    err_src(:,i) = err_src_;
end
save(save_file_name, 'err_ase', 'err_src', 'rho', 'B', 'n', 'numsim', 'K', 's_vec', 'd_vec');

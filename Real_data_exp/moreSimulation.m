% get more simulation results 
% get the table that says this is a good paper!
addpath '~/Dropbox/Research/General_code';
%nBlock = 2; 
rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
err = zeros(591, 3);
n_vec = 10:600;
parfor j = 1:length(n_vec)
    n=n_vec(j);
    [A, tau, ~] = sbm(rho, B, n);
    embed = svdembed(2, A);
    temp = zeros(1,3);
    [temp(1), ~] = errfun(tau, A, 5); 
    [temp(2), ~] = errsvdknn(tau, embed, 1);
    [temp(3), ~] = errsvdlda(tau, embed);
    err(j,:) = temp;
    
end

plot(n_vec, err(:,1), '-*k',n_vec, err(:,2), '-ok',n_vec, err(:,3),...
    '-^k', 'LineWidth', 2, 'MarkerSize',8)
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest');
xlabel('Occlusion rate','fontsize',14)
ylabel('Error rate','fontsize',14)
hold on; hline = refline([0 0.5771]);
set(hline,'Color','k')

hold off;
imagesc(A)


plot(n_vec, err(:,1), '-k',n_vec, err(:,2), '--k',n_vec, err(:,3), )
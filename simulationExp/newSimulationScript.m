% new simulation script
% in this script, we change the sparsity level
% get the table that says this is a good paper!
addpath '~/Dropbox/Research/General_code';

rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assessing Robustness
% occlusion model misspecification

n = 200;
p = 0: 0.03: 1;
dim1 = n;
s_vec = 2:5;
d_vec = [ 2 4];
err = zeros(length(p),length(s_vec));
errLDA = zeros(length(p), 2);
errNN = zeros(length(p), 2);
nsim = 1000;

parfor k= 1:length(p)     
    errTmp = zeros(nsim, length(s_vec));
    errTmpLDA = zeros(nsim, 2);
    errTmpNN = zeros(nsim, 2);
    for nexp = 1:nsim       
        % occlusion procedure
        [A, tau, ~] = sbm(rho, B, n);
        inc = round(dim1*p(k));
        i = randsample(1:dim1, inc);
        patch = A; 
        patch(i,i) = 0;
        % embed to different dimensions
        for j = 1:length(d_vec)
            d = d_vec(j);
            embed = adjSpecEmbed(patch, d);
            errTmpLDA(nexp, j) = errsvdlda(tau, embed);
            errTmpNN(nexp, j) = errsvdknn(tau, embed, 1);
        end
        % changing sparsity level
        for j = 1:length(s_vec)        
            s = s_vec(j);
            errTmp(nexp, j) = errfun(tau, patch, s);   
        end
    end
     err(k,:) = mean(errTmp);   
     errLDA(k,:) = mean(errTmpLDA);
     errNN(k,:) = mean(errTmpNN);
end

figure; plot(p(1:33), err(1:33,1), 'bo-', 'LineWidth', 2, 'MarkerSize', 6)
hold on;
plot(p(1:33), err(1:33,2), 'bd-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p(1:33), err(1:33,3), 'bs-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p(1:33), err(1:33,4), 'b^-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p(1:33), errLDA(1:33,1), 'r^-',p(1:33), errLDA(1:33, 2), 'rd-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p(1:33), errNN(1:33,1),'kx-',p(1:33), errNN(1:33,2), 'ks-', 'LineWidth', 2, 'MarkerSize', 6)
hold off;
legend('SRC$_{s=2}$', 'SRC$_{s=3}$', 'SRC$_{s=4}$', 'SRC$_{s=5}$',...
    'LDA$\circ$ASE$_{\hat{d} = 2}$', 'LDA$\circ$ASE$_{d = 4}$',...
    'NN$\circ$ASE$_{\hat{d} = 2}$', 'NN$\circ$ASE$_{d = 4}$')
set(gca,'FontSize',18)
xlabel('Occlusion rate p_o','fontsize',32)
ylabel('Error rate','fontsize',32)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Flip

n = 200;
p = 0: 0.03: 1;
dim1 = n;
s_vec = 2:5;
d_vec = [ 2 4];
err = zeros(length(p),length(s_vec));
errLDA = zeros(length(p), 2);
errNN = zeros(length(p), 2);
nsim = 1000;

parfor k= 1:length(p)     
    errTmp = zeros(nsim, length(s_vec));
    errTmpLDA = zeros(nsim, 2);
    errTmpNN = zeros(nsim, 2);
    for nexp = 1:nsim       
        % occlusion procedure
        [A, tau, ~] = sbm(rho, B, n);
        inc = round(dim1*p(k));
        i = randsample(1:dim1, inc);
        patch = A; 
        patch(i,i)= 1-patch(i,i);
        % embed to different dimensions
        for j = 1:length(d_vec)
            d = d_vec(j);
            embed = adjSpecEmbed(patch, d);
            errTmpLDA(nexp, j) = errsvdlda(tau, embed);
            errTmpNN(nexp, j) = errsvdknn(tau, embed, 1);
        end
        % changing sparsity level
        for j = 1:length(s_vec)        
            s = s_vec(j);
            errTmp(nexp, j) = errfun(tau, patch, s);   
        end
    end
     err(k,:) = mean(errTmp);   
     errLDA(k,:) = mean(errTmpLDA);
     errNN(k,:) = mean(errTmpNN);
end

figure; %plot(p, err(:,1), 'bo-', 'LineWidth', 2, 'MarkerSize', 6)
hold on;
plot(p, err(:,2), 'bd-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p, err(:,3), 'bs-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p, err(:,4), 'b^-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p, errLDA(:,1), 'r^-',p, errLDA(:, 2), 'rd-', 'LineWidth', 2, 'MarkerSize', 6)
plot(p, errNN(:,1),'kx-',p, errNN(:,2), 'ks-', 'LineWidth', 2, 'MarkerSize', 6)
hold off;
legend('SRC$_{s=2}$', 'SRC$_{s=3}$', 'SRC$_{s=4}$', 'SRC$_{s=5}$',...
    'LDA$\circ$ASE$_{\hat{d} = 2}$', 'LDA$\circ$ASE$_{d = 4}$',...
    'NN$\circ$ASE$_{\hat{d} = 2}$', 'NN$\circ$ASE$_{d = 4}$')
set(gca,'FontSize',18)
xlabel('Linkage reversion rate p_l','fontsize',32)
ylabel('Error rate','fontsize',32)


























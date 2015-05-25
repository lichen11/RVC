% get more simulation results 
% get the table that says this is a good paper!
addpath '~/Dropbox/Research/General_code';
%nBlock = 2; 
rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
n_vec = 10:100;
err = zeros(length(n_vec), 3);
parfor j = 1:length(n_vec)
    n = n_vec(j);
    [A, tau, ~] = sbm(rho, B, n);
    %A = diagAug(A);
    embed = svdembed(2, A);
    temp = zeros(1,3);
    [temp(1), ~] = errfun(tau, A, 5); 
    [temp(2), ~] = errsvdknn(tau, embed, 1);
    [temp(3), ~] = errsvdlda(tau, embed);
    err(j,:) = temp;
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
n_vec = 30:10:140;
nexp = 1750;
err = zeros(length(n_vec), 3, nexp);
for nsim = 1:nexp
    for j = 1:length(n_vec)
        n = n_vec(j);     
        [A, tau, ~] = sbm(rho, B, n);
         A = diagAug(A);
        embed = svdembed(2, A);
        temp = zeros(1,3);
        [temp(1), ~] = errfun(tau, A, 5); 
        [temp(2), ~] = errsvdknn(tau, embed, 1);
        [temp(3), ~] = errsvdlda(tau, embed);

        err(j,:, nsim) = temp;
    end
end

srcErr = err(:, 1, :);
nnErr  = err(:, 2, :);
ldaErr = err(:, 3, :);


figure; errorbar(mean(srcErr,3), std(srcErr,0,3)/sqrt(nexp), 'LineWidth', 2)
hold on;
errorbar(mean(nnErr,3), std(nnErr,0,3)/sqrt(nexp), 'black', 'LineWidth', 2)
errorbar(mean(ldaErr,3), std(ldaErr,0,3)/sqrt(nexp), 'red', 'LineWidth', 2)
hold off;
xlim([0.5, 12])
set(gca,'XTickLabel',num2str(((n_vec)')))

plot(n_vec, smooth(mean(srcErr,3)), '-bo', n_vec, smooth(mean(nnErr,3)), '-kx',n_vec, smooth(mean(ldaErr, 3)), '-r^',...,
    'LineWidth', 2, 'MarkerSize', 7)
set(gca,'FontSize',14)
xlim([n_vec(1) 120])
legend('SRC','NN o ASE','LDA o ASE','Location','NorthEast');
xlabel('Number of Vertices','FontSize',14)
ylabel('Error rate','fontsize',14)



%% Plot 



plot(n_vec(1:80), smooth(err(1:80,1)), '-bo',n_vec(1:80), smooth(err(1:80,2)), '-kx',n_vec(1:80), smooth(err(1:80,3)), '-r^',...,
    'LineWidth', 2, 'MarkerSize', 6)
plot(n_vec, err(:,1), 'bo',n_vec, err(:,2), 'kx',n_vec, err(:,3), 'r^',...,
    'LineWidth', 2, 'MarkerSize', 8)
set(gca,'FontSize',14)
legend('SRC','NN o ASE','LDA o ASE','Location','NorthEast');
xlabel('Number of Vertices','FontSize',14)
ylabel('Error rate','fontsize',14)


plot(n_vec(1:250), err(1:250,1), 'bo',n_vec(1:250), err(1:250,2), 'kx',n_vec(1:250), err(1:250,3), 'r^',...,
    'LineWidth', 2, 'MarkerSize', 8)
set(gca,'FontSize',18)
legend('SRC','1NN o ASE_d_=_2','LDA o ASE_d_=_2','Location','NorthEast');
xlabel('Number of Vertices','FontSize',32)
ylabel('Error rate','fontsize',32)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% now we vary B 

rho = [0.4 0.6];
alpha_vec = 0.1:0.05:0.91;
beta_vec = 0.1:0.05:0.91;
n = 250;
errVaryingAlphaBeta = zeros(17, 17, 3);
for j = 1:length(alpha_vec)
    for i = 1:length(beta_vec)
        alpha = alpha_vec(j);
        beta = beta_vec(i);
        B = [alpha beta; beta alpha];
        [A, tau, ~] = sbm(rho, B, n);
        embed = svdembed(2, A);
        temp = zeros(1,3);
        [temp(1), ~] = errfun(tau, A, 5); 
        [temp(2), ~] = errsvdknn(tau, embed, 1);
        [temp(3), ~] = errsvdlda(tau, embed);
        errVaryingAlphaBeta(i,j,:) = temp;
    end
end

plot(alpha_vec, errVaryingB(:,1), 'bo-',alpha_vec, errVaryingB(:,2), 'kx-',...
    alpha_vec, errVaryingB(:,3), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 8)
set(gca,'FontSize',18)
legend('SRC','NN o SVD','LDA o SVD','Location','NorthEast');
xlabel('Number of Vertices','FontSize',32)
ylabel('Error rate','fontsize',32)

% plot surface
surf(alpha_vec, beta_vec, errVaryingAlphaBeta(:,:,1))
xlabel('Alpha', 'FontSize', 20)
ylabel('Beta', 'FontSize', 20)
zlabel('Error rate', 'FontSize', 20)

% plot heatmap
imagesc(errVaryingAlphaBeta(:,:,1))
set(gca, 'YTick', alpha_vec(1:2:end))
imagesc(errVaryingAlphaBeta(:,:,1))
set(gca,'XTickLabel',num2str((alpha_vec(1:2:end)')))
set(gca,'YTickLabel',num2str((beta_vec(1:2:end)')))
xlabel('Alpha', 'FontSize', 20)
ylabel('Beta', 'FontSize', 20)

%%

b = 0.3;
rho = [0.4 0.6];
a_vec = 0.15:0.03:0.5;
n_vec = 100:50:300;
nsim = 30;
errVaryingNverAndAlphaBetaDifference = zeros(length(a_vec), 3, ...
    length(n_vec), nsim);
for nexp = 1 : nsim
    for k = 1 : length(n_vec)
        for j = 1:length(a_vec)
    
            a = a_vec(j);
            n = n_vec(k);
            B = [a+b b; b b+a];
            [A, tau, ~] = sbm(rho, B, n);
            %embed = svdembed(2, A);
            embed = adjSpecEmbed(diagAug(A), 2);
            temp = zeros(1,3);
            [temp(1), ~] = errfun(tau, A, 5); 
            [temp(2), ~] = errsvdknn(tau, embed, 1);
            [temp(3), ~] = errsvdlda(tau, embed);
            errVaryingNverAndAlphaBetaDifference(j,:, k, nexp) = temp;
        end
    end
end

srcPerform = reshape(errVaryingNverAndAlphaBetaDifference(:,1,:, :), length(a_vec), length(n_vec), nsim);
nnPerform = reshape(errVaryingNverAndAlphaBetaDifference(:,2,:, :), length(a_vec), length(n_vec), nsim);
ldaPerform = reshape(errVaryingNverAndAlphaBetaDifference(:,3,:, :), length(a_vec), length(n_vec), nsim);
srcNNDiff = srcPerform - nnPerform;
srcLDADiff = srcPerform - ldaPerform;

set(gca,'fontsize', 35) 
clims = [min(min(min(mean(srcLDADiff, 3)), min(mean(srcNNDiff, 3)))) max(max(max(mean(srcLDADiff, 3)), max(mean(srcNNDiff, 3))))];
imagesc(n_vec, a_vec, mean(srcNNDiff, 3), clims);
figure; set(gca,'fontsize', 35) ; imagesc(n_vec, a_vec, mean(srcLDADiff, 3), clims);

surf(n_vec, a_vec, mean(srcPerform, 3));
surf(srcNNDiff);
surf(srcLDADiff);
imagesc(srcLDADiff);
imagesc( n_vec, a_vec(6:length(a_vec)),srcNNDiff(6:length(a_vec),:))
imagesc(n_vec, a_vec(6:length(a_vec)), srcLDADiff(6:length(a_vec),:))

set(gca,'XTickLabel',num2str((n_vec')))
set(gca,'YTickLabel',num2str((a_vec')))
xlabel('Number of Vertices', 'FontSize', 35)
ylabel('Block Dissimilarity', 'FontSize', 35)
zlabel('Error Rate', 'FontSize', 35)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now we examine changing what happens at nVec = 100
% Is it noise or significant?
b = 0.3;
rho = [0.4 0.6];
a_vec = 0.1:0.03:0.18; % tried 0:0.02:0.18
%n_vec = 100:50:300;
nsim = 2000;
 n = 40;
errVaryingNverAndAlphaBetaDifference = zeros(length(a_vec), 3, nsim);
for k = 1 : nsim
    for j = 1:length(a_vec)    
        a = a_vec(j);
       
        B = [a+b b; b b+a];
        [A, tau, ~] = sbm(rho, B, n);
        %embed = svdembed(2, A);
        embed = adjSpecEmbed(patch, 2);
        temp = zeros(1,3);
        [temp(1), ~] = errfun(tau, A, 5); 
        [temp(2), ~] = errsvdknn(tau, embed, 1);
        [temp(3), ~] = errsvdlda(tau, embed);
        errVaryingNverAndAlphaBetaDifference(j,:,k) = temp;
    end
end
errorbar(a_vec, mean(srcPerform,2), std(srcPerform,0,2)/sqrt(nsim), 'LineWidth', 2)
hold on;
errorbar(a_vec, mean(nnPerform,2), std(nnPerform,0,2)/sqrt(nsim), 'black', 'LineWidth', 2)
errorbar(a_vec, mean(ldaPerform,2), std(ldaPerform,0,2)/sqrt(nsim), 'red', 'LineWidth', 2)
hold off;
set(gca,'XTickLabel',num2str(((a_vec)')))

srcPerform = reshape(errVaryingNverAndAlphaBetaDifference(:,1,:), length(a_vec), nsim);
nnPerform = reshape(errVaryingNverAndAlphaBetaDifference(:,2,:), length(a_vec), nsim);
ldaPerform = reshape(errVaryingNverAndAlphaBetaDifference(:,3,:), length(a_vec), nsim);
srcNNDiff = srcPerform - nnPerform;
srcLDADiff = srcPerform - ldaPerform;

pval = zeros(length(a_vec), 2);
for i = 1:length(a_vec)
    pval(i, 1) = ranksum(srcPerform(i,:), nnPerform(i,:));
    pval(i, 2) = ranksum(srcPerform(i,:), ldaPerform(i,:));
end

pval = zeros(length(a_vec), 2);
for i = 1:length(a_vec)
    [~, pval(i, 1)] = ttest(srcPerform(i,:), nnPerform(i,:));
    [~, pval(i, 2)] = ttest(srcPerform(i,:), ldaPerform(i,:));
end

errorbar(mean(srcNNDiff,2), std(srcNNDiff,0,2), 'LineWidth', 2)
hold on;
errorbar(mean(srcLDADiff,2), std(srcLDADiff,0,2), 'black', 'LineWidth', 2)
errorbar(mean(ldaPerform,2), std(ldaPerform,0,2), 'red', 'LineWidth', 2)
hold off;
xlim([0.5, 12])
set(gca,'XTickLabel',num2str(((a_vec)')))

plot(a_vec, ldaPerform);
imagesc(srcNNDiff);
imagesc(srcLDADiff);
imagesc(srcLDADiff);
imagesc( n_vec, a_vec(6:length(a_vec)),srcNNDiff(6:length(a_vec),:))
imagesc(n_vec, a_vec(6:length(a_vec)), srcLDADiff(6:length(a_vec),:))
set(gca,'fontsize', 35) 
set(gca,'XTickLabel',num2str((n_vec')))
set(gca,'YTickLabel',num2str((a_vec')))
xlabel('Number of Vertices', 'FontSize', 35)
ylabel('Block Dissimilarity', 'FontSize', 35)
zlabel('Error Rate', 'FontSize', 35)


srcPerform = srcPerform(1:10,:);
nnPerform = nnPerform(1:10, :);
ldaPerform = ldaPerform(1:10, :);
srcNNDiff = srcPerform - nnPerform;
srcLDADiff = srcPerform - ldaPerform;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dissimilarity and sparsity and nvertices
n= 450;
b= 0.4;
rho = [0.4 0.6];
a_vec = 0:0.05:0.5;
s_vec = 3:1:15;
%n_vec = 100:50:450;
errVaryingDissimAndSparsity = zeros(length(a_vec), length(s_vec));
for k = 1 : length(s_vec)
    for j = 1:length(a_vec)
        s = s_vec(k);
        a = a_vec(j);
        %n = n_vec(k);
        
        B = [a+b b; b b+a];
        [A, tau, ~] = sbm(rho, B, n);
        %embed = svdembed(2, A);
        temp = zeros(1,1);
        [temp(1), ~] = errfun(tau, A, s); 
        %[temp(2), ~] = errsvdknn(tau, embed, 1);
        %[temp(3), ~] = errsvdlda(tau, embed);
        errVaryingDissimAndSparsity(j, k) = temp;
    end
end

surf(errVaryingDissimAndSparsity)
set(gca,'YTickLabel',num2str((a_vec')))
set(gca,'XTickLabel',num2str((s_vec(1:1.7:end)')))
ylabel('Block Dissimilarity', 'FontSize', 20)
xlabel('Sparsity Level', 'FontSize', 20)
zlabel('Error Rate', 'FontSize', 20)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varying dissimilarity, n, sparsity
b= 0.4;
rho = [0.4 0.6];
a_vec = 0:0.05:0.5;
s_vec = 3:1:15;
n_vec = 100:50:450;
errVaryingDissimAndSparsityAndNvertices = zeros(length(a_vec), length(s_vec), length(n_vec));
for i = 1 : length(n_vec)
    for k = 1 : length(s_vec)
        for j = 1:length(a_vec)
        n = n_vec(i);    
        s = s_vec(k);
        a = a_vec(j);
       
        
        B = [a+b b; b b+a];
        [A, tau, ~] = sbm(rho, B, n);
        %embed = svdembed(2, A);
        temp = zeros(1,1);
        [temp(1), ~] = errfun(tau, A, s); 
        %[temp(2), ~] = errsvdknn(tau, embed, 1);
        %[temp(3), ~] = errsvdlda(tau, embed);
        errVaryingDissimAndSparsityAndNvertices(j, k, i) = temp;
        end
        
    end
end

for i = 1:8
    subplot(2,4,i, 'align');
    
    imagesc(errVaryingDissimAndSparsityAndNvertices(:,:,i))
    set(gca,'YTickLabel',num2str((s_vec(1:1.2:end)')))
    set(gca,'XTickLabel',num2str((a_vec(2:1.6:end)')))
    xlabel('Block Dissimilarity', 'FontSize', 20)
    ylabel('Sparsity Level', 'FontSize', 20)
    title(sprintf('n = %d', n_vec(i)), 'FontSize', 20)

end
h = colorbar;
set(h, 'Position', [.9314 .11 .0381 .8150])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assessing Robustness
% occlusion model misspecification
rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
%n = 200;
n = 50;
p = 0: 0.03: 0.8;
%seed = RandStream('mcg16807','Seed',8 );
%RandStream.setGlobalStream(seed)
nsim= 50;
occlusionErrTensor = zeros(length(p), 3, nsim);
occlusionErrTensor4D = zeros(length(p), 3, nsim);
parfor i = 1:nsim
    [A, tau, ~] = sbm(rho, B, n);
    occlusionErr = occlusion(p, A, tau, 2, 5);
    occlusionErr4D = occlusion(p, A, tau, 4, 6);
    occlusionErrTensor(:,:,i) = occlusionErr;
    occlusionErrTensor4D(:,:,i) = occlusionErr4D;
end



srcOccErr = occlusionErrTensor(:, 1, :);
nnOccErr  = occlusionErrTensor(:, 2, :);
ldaOccErr = occlusionErrTensor(:, 3, :);
srcOccErr = reshape(srcOccErr, length(p), nsim);
nnOccErr = reshape(nnOccErr, length(p), nsim);
ldaOccErr = reshape(ldaOccErr, length(p), nsim);

srcOccErr4D = occlusionErrTensor4D(:, 1, :);
nnOccErr4D  = occlusionErrTensor4D(:, 2, :);
ldaOccErr4D = occlusionErrTensor4D(:, 3, :);
srcOccErr4D = reshape(srcOccErr4D, length(p), nsim);
nnOccErr4D = reshape(nnOccErr4D, length(p), nsim);
ldaOccErr4D = reshape(ldaOccErr4D, length(p), nsim);

figure; errorbar(p, mean(srcOccErr4D,2), std(srcOccErr4D,0,2)/sqrt(nsim), 'LineWidth', 2)
hold on;
errorbar(p, mean(nnOccErr4D,2), std(nnOccErr4D,0,2)/sqrt(nsim), 'black', 'LineWidth', 2)
errorbar(p, mean(ldaOccErr4D,2), std(ldaOccErr4D,0,2)/sqrt(nsim), 'red', 'LineWidth', 2)
hold off;
xlim([0, 1])
%set(gca,'XTickLabel',num2str(((10:10:120)')))


plot(p, mean(srcOccErr,2), 'bo-', p, mean(nnOccErr,2), 'kx-', p, mean(ldaOccErr,2), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 6)
hold on;
plot(p, smooth(mean(nnOccErr4D,2)), 'ks-.',...
    p, smooth(mean(ldaOccErr4D,2)), 'rd-.',...,
    'LineWidth', 2, 'MarkerSize', 6)
hold off;
legend('SRC','NN o ASE_2_d','LDA o ASE_2_d', 'NN o ASE_4_d','LDA o ASE_4_d','Location','NorthWest');
figure; plot(p, smooth(mean(srcOccErr,2)), 'bo-', p, smooth(mean(nnOccErr,2)), 'kx-',...
    p, smooth(mean(ldaOccErr,2)), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 6)
figure; plot(p, smooth(mean(srcOccErr4D,2)), 'bo-', p, smooth(mean(nnOccErr4D,2)), 'kx-',...
    p, smooth(mean(ldaOccErr4D,2)), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 6)
legend('SRC','NN o ASE','LDA o ASE','Location','NorthWest');

set(gca,'FontSize',18)
xlabel('Occlusion rate p_o','fontsize',32)
ylabel('Error rate','fontsize',32)




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = [0.4 0.6];
B = [0.7 0.32; 0.32 0.75]; 
n = 200;
p = 0: 0.03: 1;
%seed = RandStream('mcg16807','Seed',8 );
%RandStream.setGlobalStream(seed)
nsim = 100;
flipErrTensor = zeros(length(p), 3, nsim);
flipErrTensor4D = zeros(length(p), 3, nsim);
parfor nexp = 1 : nsim
    [A, tau, ~] = sbm(rho, B, n);
    flipErr = flipping(p, A, tau, 2, 5);
    flipErrTensor(:, :, nexp) = flipErr;
    flipErr4D = flipping(p, A, tau, 4, 6);
    flipErrTensor4D(:, :, nexp) = flipErr4D;
end




srcFlipErr = flipErrTensor(:,1,:);
srcFlipErr = reshape(srcFlipErr, length(p),nsim);
%srcFlipErr = mean(srcFlipErr, 2);

nnFlipErr = flipErrTensor(:, 2,:);
nnFlipErr = reshape(nnFlipErr, length(p),nsim);
%nnFlipErr = mean(nnFlipErr, 2);

ldaFlipErr = flipErrTensor(:,3, :);
ldaFlipErr = reshape(ldaFlipErr, length(p),nsim);
%ldaFlipErr = mean(ldaFlipErr, 2);

figure; errorbar(p, mean(srcFlipErr,2), std(srcFlipErr,0,2)/sqrt(nsim), 'LineWidth', 2)
hold on;
errorbar(p, mean(nnFlipErr,2), std(nnFlipErr,0,2)/sqrt(nsim), 'black', 'LineWidth', 2)
errorbar(p, mean(ldaFlipErr,2), std(ldaFlipErr,0,2)/sqrt(nsim), 'red', 'LineWidth', 2)
hold off;
xlim([0, 1])
%set(gca,'XTickLabel',num2str(((10:10:120)')))


plot(p, mean(srcFlipErr,2), 'bo-', p, mean(nnFlipErr,2), 'kx-', p, mean(ldaFlipErr,2), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 6)

figure;
plot(p, smooth(mean(srcFlipErr,2)), 'bo-', p, smooth(mean(nnFlipErr,2)), 'kx-', p, smooth(mean(ldaFlipErr,2)), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 7)
set(gca,'FontSize',18)
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest');
xlabel('Linkage reversion rate p_l','fontsize',32)
ylabel('Error rate','fontsize',32)



srcFlipErr4D = flipErrTensor4D(:,1,:);
srcFlipErr4D = reshape(srcFlipErr4D, length(p),nsim);
%srcFlipErr = mean(srcFlipErr, 2);

nnFlipErr4D = flipErrTensor4D(:, 2,:);
nnFlipErr4D = reshape(nnFlipErr4D, length(p),nsim);
%nnFlipErr = mean(nnFlipErr, 2);

ldaFlipErr4D = flipErrTensor4D(:,3, :);
ldaFlipErr4D = reshape(ldaFlipErr4D, length(p),nsim);

figure; errorbar(p, mean(srcFlipErr4D,2), std(srcFlipErr4D,0,2)/sqrt(nsim), 'LineWidth', 2)
hold on;
errorbar(p, mean(nnFlipErr4D,2), std(nnFlipErr4D,0,2)/sqrt(nsim), 'black', 'LineWidth', 2)
errorbar(p, mean(ldaFlipErr4D,2), std(ldaFlipErr4D,0,2)/sqrt(nsim), 'red', 'LineWidth', 2)
hold off;
xlim([0, 1])
%set(gca,'XTickLabel',num2str(((10:10:120)')))


plot(p, mean(srcFlipErr,2), 'bo-', p, mean(nnFlipErr,2), 'kx-', p, mean(ldaFlipErr,2), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 6)

figure;
plot(p, smooth(mean(srcFlipErr,2)), 'bo-', p, smooth(mean(nnFlipErr,2)), 'kx-', p, smooth(mean(ldaFlipErr,2)), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 7)
hold on;
plot(p, smooth(mean(nnFlipErr4D,2)), 'ks-.', p, smooth(mean(ldaFlipErr4D,2)), 'rd-.',...,
    'LineWidth', 2, 'MarkerSize', 7)
h = legend('SRC','1NN$\circ$ASE$_{\hat{d}=2}$','LDA$\circ$ASE$_{\hat{d}=2}$', '1NN$\circ$ASE$_{d=4}$',...
    'LDA$\circ$ASE$_{d=4}$','Location','NorthWest');
set(h,'Interpreter','latex')
hold off;

plot(p, smooth(mean(srcFlipErr4D,2)), 'go-', p, smooth(mean(nnFlipErr4D,2)), 'kx-', p, smooth(mean(ldaFlipErr4D,2)), 'r^-',...,
    'LineWidth', 2, 'MarkerSize', 7)


%[A, tau, ~] = sbm(rho, B, n);
for k= 1:length(p)
  inc = round(dim1*p(k));
  %seed = rng;
  i = randsample(1:dim1, inc);

  patch1 = adj; 
  patch1(i,i) = 0;  %occlusion
  %patch1 = diagAug(patch1);
  [~, D1] = eigs(patch1, 4);
  
  patch2 = adj; 
  patch2(i,i)= 1-patch2(i,i); 
  patch2 = diagAug(patch2);
  [~, D2] = eigs(patch2, 4); 
end

p = 0:0.01:1;
dim1 = n;
[adj, tau, ~] = sbm(rho, B, n);


set(gca,'FontSize',18)
for k = 1:length(p)
  inc = round(dim1*p(k));
  %seed = rng;
  i = randsample(1:dim1, inc);
  patch1 = adj; 
  patch1(i,i)= 0; 
  %patch1 = diagAug(patch1);
  [~, D1] = eigs(patch1, 22);
  D1 = diag(D1);
  group = (D1>0);
  tmp = [group D1 abs(D1)];
  tmp = -sortrows(-tmp, 3);
  
  gscatter(1:22, tmp(:,3), tmp(:,1), ['r', 'g'], [], 28)
  %Dposi = D1(D1>0);
  %Dneg = D1(D1<0);
  %plot(abs(Dposi), 'ro', 'go', 'LineWidth', 2,'MarkerSize', 8);
  %hold on;
  %plot(abs(Dneg), 'go', 'LineWidth', 2,'MarkerSize', 8);
  %hold off;
  title(sprintf('p=%g', p(k)))
  set(gca,'FontSize',18)
  legend('Negative','Positive','Location','NorthEast');
  xlim([0 22])   
  ylim([0 120])
  set(gca,'XTickLabel',[])
  M(k) = getframe(gcf);
  pause(0.01);  
end
movie2avi(M, 'eVal22OccChange.avi');

set(gca,'FontSize',18)
for k = 1:length(p)
  inc = round(dim1*p(k));
  %seed = rng;
  i = randsample(1:dim1, inc);

  patch1 = adj; 
  patch1(i,i)= 1-patch1(i,i); 
  %patch1 = diagAug(patch1);
  [~, D1] = eigs(patch1, 4);
  plot(abs(diag(D1)), 'o',  'LineWidth', 2,'MarkerSize', 8)
  title(sprintf('p=%g', p(k)))
  xlim([0 5])
   
  ylim([0 120])
  set(gca,'XTickLabel',[])
  M(k) = getframe(gcf);
  pause(0.02); 
 
end





p = 0.74;
d_vec = 1:125;
nsim = 150;
errorMat = zeros(length(d_vec), 3);
for d = d_vec
    occlusionErrMat = zeros(nsim, 3);
    parfor i = 1:nsim
        [A, tau, ~] = sbm(rho, B, n);
        occlusionErr = occlusion(p, A, tau, d, d);  
        occlusionErrMat(i,:) = occlusionErr;
    end
    errorMat(d,:) = mean(occlusionErrMat,1);
end

save('changeDandSresultFrom1To125.mat', 'd_vec', 'p', 'nsim', 'errorMat');

figure; plot(1:(d-1), errorMat(1:(d-1),1), 'bo-', 1:(d-1), errorMat(1:(d-1),2), 'kx-',1:(d-1), errorMat(1:(d-1),3), 'r^-' ,'LineWidth', 2, 'MarkerSize', 6);
set(gca,'FontSize',18)
legend('SRC','1NN o ASE','LDA o ASE','Location','NorthEast');
xlabel('d (or s)','FontSize',32)
ylabel('Error rate','fontsize',32)

nsim = 150;
errLDA = zeros(length(103:130),1);
errNN = zeros(length(103:130), 1);
dim1 = n;
for d = 103:130
    errLDAtmp = zeros(nsim,1);
    errNNtmp = zeros(nsim, 1);
    parfor j = 1:nsim
        [A, tau, ~] = sbm(rho, B, n);
        
        inc = round(dim1*p);
        i = randsample(1:dim1, inc);

        patch = A; 
        patch(i,i) = 0;  %occlusion

        embed = adjSpecEmbed(patch, d);
        [errLDAtmp(j), ~] = errsvdlda(tau, embed);
        [errNNtmp(j), ~] = errsvdknn(tau, embed, 1);
    end
    errLDA(d-102) = mean(errLDAtmp);
    errNN(d-102) = mean(errNNtmp);
end

hold on; 
plot(103:127, errNN(1:25),'kx-', 103:127, errLDA(1:25), 'r^-' ,'LineWidth', 2, 'MarkerSize', 6);
set(gca,'FontSize',18)



tic;
nsim = 150;
d_vec = 1:3:89;
errLDA = zeros(length(d_vec),1);
errNN = zeros(length(d_vec), 5);
errSRC = zeros(length(d_vec), 1);
%err5NN = zeros(length(d_vec), 1);

for i = 1:length(d_vec)
    d = d_vec(i);
    errLDAtmp = zeros(nsim,1);
    errNNtmp = zeros(nsim, 5);
    %errSRCtmp = zeros(nsim, 1);
    parfor j = 1:nsim
        [A, tau, ~] = sbm(rho, B, n);
        embed = adjSpecEmbed(A, d);
      %  errSRCtmp(j) = errfun(tau, A, d);
        [errLDAtmp(j), ~] = errsvdlda(tau, embed);
        
        [tmp1, ~] = errsvdknn(tau, embed, 1);
        [tmp2, ~] = errsvdknn(tau, embed, 3);

        [tmp3, ~] = errsvdknn(tau, embed, 9);
        [tmp4, ~] = errsvdknn(tau, embed, 13);
        [tmp5, ~] = errsvdknn(tau, embed, 17);
        errNNtmp(j, :) = [tmp1 tmp2 tmp3 tmp4 tmp5];
       
        
     %   [err3NNtmp(j), ~] = errsvdknn(tau, embed, 3);
     %   [err5NNtmp(j), ~] = errsvdknn(tau, embed, 5);

    end
    errLDA(i) = mean(errLDAtmp);
    errNN(i,:) = mean(errNNtmp,1);
    %errSRC(i) = mean(errSRCtmp);
    %err3NN(d) = mean(err3NNtmp);
    %err5NN(d) = mean(err5NNtmp);
end
toc;
figure; plot( d_vec, errLDA, 'r^-', d_vec, errNN, 'kx-', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(d_vec, errSRC, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
hold off;


figure; plot( d_vec, errLDA, 'r^-',  'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(d_vec, errNN(:,1), 'kx-','LineWidth', 2, 'MarkerSize', 6);
plot(d_vec, errNN(:,2), 'gx-',  'LineWidth', 2, 'MarkerSize', 6);
plot(d_vec, errNN(:,3), 'bx-',  'LineWidth', 2, 'MarkerSize', 6);
plot(d_vec, errNN(:,4), 'k^-',  'LineWidth', 2, 'MarkerSize', 6);
plot(d_vec, errNN(:,5), 'yx-',  'LineWidth', 2, 'MarkerSize', 6);
plot(d_vec, errSRC, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
hold off;

legend('LDA', '1NN','3NN','9NN','13NN','17NN')
save('LDAandNN', 'errLDA', 'd_vec', 'errNN')


%load run_SRC_on_Football_1_100
load run_SRC_on_Adjnoun_1_100
%load run_SRC_on_Blog_1_100
%load run_SRC_on_Polibook_1_100
%load run_SRC_on_celegans_1_100 
subind = 1:50;
plot(d_vec(subind), errMat(subind,1), 'bo-','LineWidth', 3, 'MarkerSize', 8);
hold on;
plot(d_vec(subind), errMat(subind,2), 'gx-','LineWidth', 3, 'MarkerSize', 8);
plot(d_vec(subind), errMat(subind,3), 'r^-','LineWidth', 3, 'MarkerSize', 8);
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([subind(1) subind(end) 0.1 0.53])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthWest');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Classification on AdjNoun Graph', 'FontSize', 20)









load resultBetaCEAgSummary
d_vec = 1:sequential;
plot(d_vec, mean(eomp, 2), 'bo-','LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(d_vec, mean(eASEKNN,2), 'kx-','LineWidth', 2, 'MarkerSize', 8);
plot(d_vec, mean(eASELDA,2),'r^-','LineWidth', 2, 'MarkerSize', 8);
hold off;
legend('SRC$_{s}$','9NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Classification on Celegans Ag Graph', 'FontSize', 20)
set(gca,'FontSize',20, 'fontWeight','bold')


% plot all stuff

d_vec = 1:sequential;
plot(d_vec, mean(eomp, 2), 'bo-','LineWidth', 2, 'MarkerSize', 7);
hold on;
plot(d_vec, mean(el1, 2), 'md-','LineWidth', 2, 'MarkerSize', 7); %el1 homotopy
plot(d_vec, mean(eknn, 2), 'cs-','LineWidth', 2, 'MarkerSize', 7) % marginal regression
plot(d_vec, mean(eASEKNN,2), 'kx-','LineWidth', 2, 'MarkerSize', 8);
plot(d_vec, mean(eASELDA,2),'r^-','LineWidth', 2, 'MarkerSize', 8);

hold off;
legend('SRC$_{s}$(OMP)', 'SRC$_{s}$(L$_1$)', 'SRC$_{s}$(Marginal Regression)', ...
    '9NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$', ...
    'Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Classification on Celegans Ag Graph', 'FontSize', 20)
set(gca,'FontSize',20, 'fontWeight','bold')


% plot enron
addpath '../RealDataExpNew/enronResults/';
load enronVaryingSparsityLevelFrom1To80
load enronVaryingEmbeddingDimFrom1To80NewLabel.mat
d_vec = 1:80;
plot(d_vec, err, 'bo-','LineWidth', 3, 'MarkerSize', 8);
hold on;
plot(d_vec, errVaryingD(:,1), 'gx-','LineWidth', 3, 'MarkerSize', 8);
plot(d_vec, errVaryingD(:,2),'r^-','LineWidth', 3, 'MarkerSize', 8);
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([0.5 80 0 0.46])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$', ...
    'Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Classification on Enron Graph', 'FontSize', 20)

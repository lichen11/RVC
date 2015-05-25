%plot no contamination: adjnoun like simulation
load simulation_adjnounlike_n_120_no_contamination_s_d
mean_err = mean(err, 3);
plot(d_vec, mean_err(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(d_vec, mean_err(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(d_vec, mean_err(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([d_vec(1) d_vec(end) 0.09 0.5])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','Northwest');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('No contamination: n = 110', 'FontSize', 20)


load simulation_adjnounlike_no_contamination_vary_n_50_350.mat
sub = 1:length(n_vec);
mean_err = mean(err,3);
plot(n_vec(sub), mean_err(sub,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(n_vec(sub), mean_err(sub,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(n_vec(sub), mean_err(sub,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([n_vec(1) n_vec(end) 0 0.45])
legend('SRC$_{5}$','1NN$\circ$ASE$_{d = 2}$','LDA$\circ$ASE$_{d=2}$','Location','NorthEast');
xlabel('Number of vertices','FontSize',20)
ylabel('Error rate','fontsize',20)
title('No contamination', 'FontSize', 20)


plot(n_vec, mean_err(:,1), 'bo-','LineWidth', 2, 'MarkerSize', 10)
hold on;
plot(n_vec, mean_err(:,2), 'kx-','LineWidth', 2, 'MarkerSize', 10)
plot(n_vec, mean_err(:,3), 'r^-','LineWidth', 2, 'MarkerSize', 10)
hold off;
axis([n_vec(1) n_vec(end) 0 0.45])
legend('SRC$_{5}$','1NN$\circ$ASE$_{d = 2}$','LDA$\circ$ASE$_{d=2}$','Location','NorthEast');
xlabel('Number of vertices','FontSize',20)
ylabel('Error rate','fontsize',20)
title('No contamination', 'FontSize', 20)
set(gca,'FontSize',20, 'fontWeight','bold')


%plot no contamination
load simulation_no_contamination
mean_err = mean(err,3);
plot(n_vec, mean_err(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(n_vec, mean_err(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(n_vec, mean_err(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([n_vec(1) n_vec(end) 0 0.25])
legend('SRC$_{5}$','1NN$\circ$ASE$_{d = 2}$','LDA$\circ$ASE$_{d=2}$','Location','NorthEast');
xlabel('Number of vertices','FontSize',20)
ylabel('Error rate','fontsize',20)
title('No contamination', 'FontSize', 20)

% error plot
std_err = 1/sqrt(numsim)*std(err, 0 ,3);
errorbar(n_vec, mean_err(:,1), std_err(:, 1),'bo-','LineWidth', 2, 'MarkerSize', 8); 
hold on;
errorbar(n_vec, mean_err(:, 2), std_err(:, 2),'kx-','LineWidth', 2, 'MarkerSize', 8 ); 
errorbar(n_vec, mean_err(:, 3), std_err(:, 3),'r^-','LineWidth', 2, 'MarkerSize', 8 ); 
hold off;
axis([n_vec(1) n_vec(end) 0 0.25])

% plot no contamination
%load simulation_n_40_no_contamination_s_d
load simulation_no_contamination_s_d
mean_err = mean(err,3);
plot(d_vec, mean_err(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(d_vec, mean_err(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(d_vec, mean_err(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
%axis([d_vec(1) d_vec(end) 0 0.25])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('No contamination: n = 110', 'FontSize', 20)
%
figure; 
errase = mean(err_ase,3);
errsrc = mean(err_src,2);
plot(s_vec, errsrc, 'bo-','LineWidth', 2, 'MarkerSize', 10)
hold on;
plot(d_vec, errase(:,1), 'kx-','LineWidth', 2, 'MarkerSize', 10)
plot(d_vec, errase(:,2), 'r^-','LineWidth', 2, 'MarkerSize', 10)
hold off;
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
axis([d_vec(1) d_vec(end) 0 0.112])
title('No contamination: n = 100', 'FontSize', 20)
set(gca,'FontSize',20, 'fontWeight','bold')

% plot occlusion
load run_occlusion_sparsity_embedding
mean_err1 = mean(occlusionErrTensor,3);
mean_err2 = mean(occlusionErrTensor4D, 3);
plot(p, mean_err1(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(p, mean_err1(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err1(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,2), 'gs-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,3), 'rd-','LineWidth', 3, 'MarkerSize', 10)

hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([p(1) p(end) 0 0.3])
legend('SRC$_5$','1NN$\circ$ASE$_{\hat{d}=2}$','LDA$\circ$ASE$_{\hat{d}=2}$','1NN$\circ$ASE$_{d=4}$','LDA$\circ$ASE$_{d=4}$','Location','NorthWest');
xlabel('Occlusion rate','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Occlusion Contamination', 'FontSize', 20)



% plot occlusion 
load occlusion_pt6_vary_s_d
mean_acc = mean(occlusionErrTensor, 3); % oops acc is misnamed, should be err
%figure;
plot(d_vec, mean_acc(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(d_vec, mean_acc(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(d_vec, mean_acc(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([d_vec(1) d_vec(end) 0 0.6])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Occlusion rate 0.6', 'FontSize', 20)

subind = 1:5;
d_vec_sub = d_vec(subind);
std_err = 1/sqrt(nsim)*std(occlusionErrTensor, 0 ,3);
errorbar(d_vec_sub, mean_err(subind, 1), std_err(subind, 1),'bo-','LineWidth', 2, 'MarkerSize', 8 ); 
hold on;
errorbar(d_vec_sub, mean_err(subind, 2), std_err(subind, 2),'kx-','LineWidth', 2, 'MarkerSize', 8 ); 
errorbar(d_vec_sub, mean_err(subind, 3), std_err(subind, 3),'r^-','LineWidth', 2, 'MarkerSize', 8 ); 
hold off;


figure;
plot(d_vec_sub, mean_acc(subind,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(d_vec_sub, mean_acc(subind,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(d_vec_sub, mean_acc(subind,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
axis([d_vec_sub(1) d_vec_sub(end) 0 0.7])
set(gca,'FontSize',20, 'fontWeight','bold')
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Occlusion rate 0.9', 'FontSize', 20)

%plot occlusion vary s k
load occlusion_pt9_vary_s_k_d4
d4_err = occlusionErrTensor;
load occlusion_pt9_vary_s_k_d2
d2_err = occlusionErrTensor;
mean_err2 = mean(d2_err,3);
mean_err = mean(d4_err, 3); % oops acc is misnamed, should be err
%figure;
plot(k_vec, mean_err(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(k_vec, mean_err(:,2), 'gs-','LineWidth', 3, 'MarkerSize', 10)
plot(k_vec, mean_err2(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
%plot(d_vec, mean_err(:,3), 'r^-','LineWidth', 2, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([k_vec(1) k_vec(end) 0 0.3])
legend('SRC$_{s}$','$k$NN$\circ$ASE$_{d = 4}$','$k$NN$\circ$ASE$_{\hat{d}=2}$','Location','NorthWest');
%xlabel({'Number of nearest neighbor','sparsity level'},'FontSize',20)
xlabel('Number of nearest neighbor/s','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Occlusion rate 0.9', 'FontSize', 20)

% reversion
load run_reversion_sparsity_embedding
mean_err1 = mean(flipErrTensor,3);
mean_err2 = mean(flipErrTensor4D, 3);
plot(p, mean_err1(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(p, mean_err1(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err1(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,2), 'gs-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,3), 'rd-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([p(1) p(end) 0 0.18])
legend('SRC$_5$','1NN$\circ$ASE$_{\hat{d}=2}$','LDA$\circ$ASE$_{\hat{d}=2}$','1NN$\circ$ASE$_{d=4}$','LDA$\circ$ASE$_{d=4}$','Location','NorthWest');
xlabel('Linkage reversion rate','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Linkage Reversion Contamination', 'FontSize', 20)

% reversion linkage
load reversion_pt7_vary_s_d
mean_err = mean(flipErrTensor, 3); % oops acc is misnamed, should be err
%figure;
plot(d_vec, mean_err(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(d_vec, mean_err(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(d_vec, mean_err(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([d_vec(1) d_vec(end) 0 0.5])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Linkage reversion rate 0.7', 'FontSize', 20)


subind = 1:5;
d_vec_sub = d_vec(subind);
plot(d_vec_sub, mean_err(subind,1), 'bo-','LineWidth', 2, 'MarkerSize', 10)
hold on;
plot(d_vec_sub, mean_err(subind,2), 'kx-','LineWidth', 2, 'MarkerSize', 10)
plot(d_vec_sub, mean_err(subind,3), 'r^-','LineWidth', 2, 'MarkerSize', 10)
hold off;
axis([d_vec_sub(1) d_vec_sub(end) 0 0.5])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Reversion rate 0.8', 'FontSize', 20)
set(gca,'FontSize',20, 'fontWeight','bold')


% plot mix contamination
load run_occlusion_reversion
mean_err1 = mean(occlusion_reversion_ErrTensor1, 3); 
mean_err2 = mean(occlusion_reversion_ErrTensor2, 3); 
plot(p_vec, mean_err1(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(p_vec, mean_err1(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(p_vec, mean_err1(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
plot(p_vec, mean_err2(:,2), 'gs-','LineWidth', 3, 'MarkerSize', 10)
plot(p_vec, mean_err2(:,3), 'rd-','LineWidth', 3, 'MarkerSize', 10)
set(gca,'FontSize',20, 'fontWeight','bold')
axis([p_vec(1) p_vec(end) 0 0.55])
legend('SRC$_5$','1NN$\circ$ASE$_{\hat{d}=2}$','LDA$\circ$ASE$_{\hat{d}=2}$','1NN$\circ$ASE$_{d=4}$','LDA$\circ$ASE$_{d=4}$','Location','NorthWest');
xlabel('Contamination rate','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Mix Contamination', 'FontSize', 20)


% plot mix contamination vary s d
load run_occlusion_reversion_pt8_vary_s_d
mean_err = mean(occlusion_reversion_ErrTensor, 3); % oops acc is misnamed, should be err
%figure;
plot(d_vec, mean_err(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(d_vec, mean_err(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(d_vec, mean_err(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([d_vec(1) d_vec(end) 0 0.55])
legend('SRC$_{s}$','1NN$\circ$ASE$_{\hat{d}}$','LDA$\circ$ASE$_{\hat{d}}$','Location','NorthEast');
xlabel('Embedding dimension/sparsity level','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Contamination rate 0.8', 'FontSize', 20)


% plot adjnoun like simulation
load run_occlusion_sparsity_embedding_adjnoun_like
mean_err1 = mean(occlusionErrTensor,3);
mean_err2 = mean(occlusionErrTensor4D, 3);
plot(p, mean_err1(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(p, mean_err1(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err1(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,2), 'gs-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,3), 'rd-','LineWidth', 3, 'MarkerSize', 10)

hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([p(1) p(end) 0.04 0.5])
legend('SRC$_5$','1NN$\circ$ASE$_{\hat{d}=2}$','LDA$\circ$ASE$_{\hat{d}=2}$','1NN$\circ$ASE$_{d=4}$','LDA$\circ$ASE$_{d=4}$','Location','NorthWest');
xlabel('Occlusion rate','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Occlusion Contamination', 'FontSize', 20)


% plot occlusion vary s and d
load occlusion

% plot reversion
load run_reversion_sparsity_embedding_adjnoun_like
mean_err1 = mean(flipErrTensor,3);
mean_err2 = mean(flipErrTensor4D, 3);
plot(p, mean_err1(:,1), 'bo-','LineWidth', 3, 'MarkerSize', 10)
hold on;
plot(p, mean_err1(:,2), 'gx-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err1(:,3), 'r^-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,2), 'gs-','LineWidth', 3, 'MarkerSize', 10)
plot(p, mean_err2(:,3), 'rd-','LineWidth', 3, 'MarkerSize', 10)
hold off;
set(gca,'FontSize',20, 'fontWeight','bold')
axis([p(1) p(end) 0.03 0.7])
legend('SRC$_5$','1NN$\circ$ASE$_{\hat{d}=2}$','LDA$\circ$ASE$_{\hat{d}=2}$','1NN$\circ$ASE$_{d=4}$','LDA$\circ$ASE$_{d=4}$','Location','NorthWest');
xlabel('Linkage reversion rate','FontSize',20)
ylabel('Error rate','fontsize',20)
title('Linkage Reversion Contamination', 'FontSize', 20)








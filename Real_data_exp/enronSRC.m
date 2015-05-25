addpath '~/Dropbox/Research/General_code';
addpath 'enron';
load('~/Dropbox/Research/Robust_est/Real_data_exp/enron/adjAndTauOverAllTimes.mat') % change config at different comps

imagesc(1-sortedAall); colormap(pink)
hold on; refline([0 114]);
line([114 ; 114],[0 ; 184])
hold off;

A = sortedAall;
Lab = sortedTau;


%At what level of sparsity?
err = zeros();
parfor s = 1: 80
     tmp = errfun(Lab, A, s);
     err(s) = tmp;
end
figure; plot(1:80, err, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sparsity level'); ylabel('Error rate'); 

embed = svdembed(2, A);
err = zeros(3,1);
[err(1), ~] = errfun(Lab, A, 5);
[err(2), ~] = errsvdknn(Lab, embed, 1);
[err(3), ~] = errsvdlda(Lab, embed);



p = 0: 0.01: 0.95;
OccEr = occlusion(p, A, Lab, 5);

plot(p, smooth(OccEr(:,1)), '-*b',p, smooth(OccEr(:,2)), '-ok',p, smooth(OccEr(:,3)),...
    '-^r', 'LineWidth', 2, 'MarkerSize',7)
plot(0:0.01:0.94, smooth(OccEr(1:95,1)), '-*b',0:0.01:0.94, smooth(OccEr(1:95,2)), '-ok',0:0.01:0.94, smooth(OccEr(1:95,3)),...
    '-^r', 'LineWidth', 2, 'MarkerSize',7)
set(gca,'FontSize',18)
legend('SRC','NN o ASE','LDA o ASE','Location','SouthWest');
xlabel('Occlusion rate','fontsize',32)
ylabel('Error rate','fontsize',32)
%hold on; hline = refline([0 0.41]);
%set(hline,'Color','k')

%hold off;
p = 0: 0.01: 1;
FlipEr = flipping(p, A, Lab, s);
figure;

figure; plot(p, smooth(FlipEr(:,1)), '-*b',p, smooth(FlipEr(:,2)), '-ok',p, smooth(FlipEr(:,3)),...
    '-^r', 'LineWidth', 2, 'MarkerSize',7)
set(gca,'FontSize',18)
legend('SRC','NN o ASE','LDA o ASE','Location','SouthWest');
xlabel('Linkage Reversion rate','fontsize',32)
ylabel('Error rate','fontsize',32)








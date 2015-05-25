addpath '~/Dropbox/Research/General_code';
addpath 'Celegans';
addpath 'Space_archaeology';
load 12ClsMat; load 13ClsMat; load 23ClsMat;
load label; load herm_connectomes;

%First check the figure with the paper 
paperLab = label + 3; %class 1 is 4, class 2 is 5, class 3 is 6;
%class 1 is class 3 in paper, class 2 is class 1 in paper, class 3 is class
%2 in paper
paperLab(paperLab == 4) = 3;
paperLab(paperLab == 5) = 2;
paperLab(paperLab == 6) = 1;

[sortedA, sortedLab] = sortAdj(full(Agap), label);
sortedA = binarA(sortedA);


OneThreeA = binarA(OneThreeA);
OneTwoA = binarA(OneTwoA);
TwoThreeA = binarA(TwoThreeA);
BinarAgap = binarA(Agap);

%A = full(BinarAgap);
%A = full(Agap);
A = full(sortedA);
Lab = sortedLab;



%At what level of sparsity?
err = zeros();
for s = 1: 250
    err(s) = errfun(Lab, A, s);
end
figure; plot(1:10, err, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sparsity level'); ylabel('Error rate'); 

embed = svdembed(3, A);
err = zeros(3,1);
[err(1), ~] = errfun(Lab, A, 5);
[err(2), ~] = errsvdknn(Lab, embed, 1);
[err(3), ~] = errsvdlda(Lab, embed);



p = 0: 0.01: 0.86;
OccEr = occlusion(p, A, Lab, 5);


plot(p, smooth(OccEr(:,1)), '-*b',p, smooth(OccEr(:,2)), '-ok',p, smooth(OccEr(:,3)),...
    '-^r', 'LineWidth', 2, 'MarkerSize',7)
legend('SRC','NN o ASE','LDA o ASE','Location','SouthWest');
set(gca,'FontSize',18)
xlabel('Occlusion rate','fontsize',32)
ylabel('Error rate','fontsize',32)

hold on; hline = refline([0 0.5771]);
set(hline,'Color','k')

hold off;
p = 0: 0.01: 1;

FlipEr = flipping(p, A, Lab, s);
figure;
plot(p, smooth(FlipEr(:,1)), '-*b',p, smooth(FlipEr(:,2)), '-ok',p, smooth(FlipEr(:,3)),...
    '-^r', 'LineWidth', 2, 'MarkerSize', 7)
set(gca,'FontSize',18)
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest');
xlabel('Linkage Reversion rate','fontsize', 32)
ylabel('Error rate','fontsize', 32)
hold on; hline = refline([0 0.5771]);
set(hline,'Color','k')

%M = load 'ERegion_ASCII_Found_Sites.xls';
clear Northing Easting Site;
Raw = [Row Column];
DisMat = dist(Raw', 'euclidean');
mean(DisMat(find(DisMat ~= 0)))





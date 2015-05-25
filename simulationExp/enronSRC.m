addpath '~/Dropbox/Research/General_code';
addpath 'enron';
addpath 'General_code';

fid = fopen('./enron/A58-184x184.txt');
A1 = fscanf(fid, '%d %d', [184, 184]);
fclose(fid);

fid = fopen('./enron/A132-184x184.txt');
A2 = fscanf(fid, '%d %d', [184, 184]);
fclose(fid);

fid = fopen('./enron/A136-184x184.txt');
A3 = fscanf(fid, '%d %d', [184, 184]);
fclose(fid);

fid = fopen('./enron/A146-184x184.txt');
A4 = fscanf(fid, '%d %d', [184, 184]);
fclose(fid);


[sortedA, sortedLab] = sortAdj(full(Agap), paperLab);
sortedA = binarA(sortedA);


OneThreeA = binarA(OneThreeA);
OneTwoA = binarA(OneTwoA);
TwoThreeA = binarA(TwoThreeA);
BinarAgap = binarA(Agap);

%A = full(BinarAgap);
%A = full(Agap);
A = full(Achem);
Lab = label;



%At what level of sparsity?
err = zeros();
for s = 1: 10
    err(s) = errfun(Lab, A, s);
end
figure; plot(1:10, err, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sparsity level'); ylabel('Error rate'); 

embed = svdembed(3, A);
err = zeros(3,1);
[err(1), ~] = errfun(Lab, A, 5);
[err(2), ~] = errsvdknn(Lab, embed, 1);
[err(3), ~] = errsvdlda(Lab, embed);



p = 0: 0.05: 0.75;
OccEr = occlusion(p, A, Lab, 5);
%x = mean(err_tensor2, 3);
%tmp = 0.1:0.05:0.95;
%plot(tmp, x(:,1), '-*k',tmp, x(:,2), '-ok',tmp, x(:,3),...
%    '-^k', 'LineWidth', 2, 'MarkerSize',8)
plot(p, OccEr(:,1), '-*k',p, OccEr(:,2), '-ok',p, OccEr(:,3),...
    '-^k', 'LineWidth', 2, 'MarkerSize',8)
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest');
xlabel('Occlusion rate','fontsize',14)
ylabel('Error rate','fontsize',14)
hold on; hline = refline([0 0.5771]);
set(hline,'Color','k')

hold off;
p = 0: 0.05: 1;

FlipEr = flipping(p, A, Lab, s);
figure;
plot(p, FlipEr(:,1), '-*k',p, FlipEr(:,2), '-ok',p, FlipEr(:,3),...
    '-^k', 'LineWidth', 2, 'MarkerSize',8)
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest');
xlabel('Linkage Reversion rate','fontsize',14)
ylabel('Error rate','fontsize',14)
hold on; hline = refline([0 0.5771]);
set(hline,'Color','k')

%M = load 'ERegion_ASCII_Found_Sites.xls';
clear Northing Easting Site;
Raw = [Row Column];
DisMat = dist(Raw', 'euclidean');
mean(DisMat(find(DisMat ~= 0)))





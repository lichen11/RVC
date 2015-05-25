% Eigenface
% input: test index: randsample(1:2414,500)
%        d: dimension
%        imgdat: image data
%        facecls: face labels
% output: acc: accuracy 
%         pred: predicted labels
%         COEFF: PCA coeff
%         latent: e-values
%         d: the dimension for 95% variation capture
function [pred, acc, rep_test, rep_train, d] = eigenface(test, train, d, test_data, train_data, facecls)

%[~, nimg]=size(imgdat);
%nimg=1:nimg;
%train = nimg(setdiff(1:length(nimg),test));
[~, ntrain] = size(train_data);
ave = 1/ntrain*sum(train_data,2);

%PCA
tic;
[COEFF,SCORE,latent] = princomp(train_data');
toc;

while sum(latent(1:d))/sum(latent)<0.95
    d = d+10;        
end

sprintf('actual variation percentage: %d',sum(latent(1:d))/sum(latent))

rep_test = transpose(COEFF(:,1:d))*(bsxfun(@minus,...
    test_data,ave));

rep_train = SCORE(:,1:d)';

pred = knnclassify(rep_test', rep_train', facecls(train), 10);
acc = length(find((pred - facecls(test))==0))/length(test);





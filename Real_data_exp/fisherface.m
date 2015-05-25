function [class, acc] = fisherface(test, train, rep_test, rep_train, facecls)
%Input: test, train indices
%       rep_train: training imgs represented in the eigenface domain
%       rep_test: test images represented in the eigenface domain
%       facecls: face labels
%Output: class: predicted labels
%        acc: accuracy

class = classify(rep_test', rep_train', facecls(train),'linear');
acc = length(find((class - facecls(test))==0))/length(test);
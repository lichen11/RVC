function [pred, acc] = sparse_rep(test, train, test_data, train_data, facecls)

% Input: test index, train index 
%        test_data, train data
%        facecls, face labels

% Output: pred, predicted value
%         acc, accurary

%L2-normalize the columns of the training dictionary
facecls = double(facecls);
Phi = bsxfun(@times,train_data, 1./sqrt(sum(train_data.^2, 1))); 

Y = test_data; %test images

X = zeros(length(train), length(test)); %recovered test images


nlabs =length(unique(facecls));
pred = zeros(length(test),1);
D_lab = [facecls(train)'; Phi];
%Solve for L1 minimization via OMP
for i = 1:length(test)
    X(:,i) = OMP(Y(:,i), Phi, 8);
    r = zeros(nlabs,1);  %residual to each class
    for j = 1: nlabs
        index = find(D_lab(1,:) == j); 
        A = D_lab(2:size(D_lab,1), index);
        %A = Phi(1:19,index);
        r(j) = norm(A*X(index,i)-Y(:,i),2);
    end
    [~,pred(i)] = min(r);
    
end
acc = length(find((facecls(test) - pred)==0))/length(test);
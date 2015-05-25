function [OMP_pred, selected_vertices, corresponding_cls, coefficients, r, recoveredSignal] = srcRep(test, test_data, train_data, tau, s)

% Input: test index, test vertex, 
%        train_data, label(tau), s
%        

% Output: OMP predicted value, vertices that are seletected, coefficients
% weighing the vertices

%L2-normalize the columns of the training dictionary
tau = double(tau);
nimg = length(tau);               %number of vertices
nimg = 1:nimg;
train = nimg(setdiff(1:length(nimg),test)); %training index set

%L2 normalize the training data
Phi = bsxfun(@times,train_data, 1./sqrt(sum(train_data.^2, 1))); 
Phi(isnan(Phi)) = 0 ;

%Phi = P(train,  train);
%Suppose no L2 normalization
%Phi = bsxfun(@times,train_data, 1); 
%sort the training data by the labels
D_lab = [tau(train)'; Phi];
[~, d2] = sort(D_lab(1,:));
D_lab = D_lab(:,d2);
Phi = D_lab(2:size(D_lab,1),:);

Y = test_data; 
X = zeros(length(train), 1); %recovered test vertex
%X2 = zeros(length(train),1);

nlabs =length(unique(tau));
labcls = unique(tau);
%pred = zeros(length(test),1);

%Solve for L1 minimization via OMP
%for i = 1:length(test)
    X(:,1) = OMP(Y(:,1), Phi, s);
    selected_vertices = find(X ~= 0);
    corresponding_cls = tau(find(X ~= 0));
    coefficients = X(X ~= 0);
    recoveredSignal = Phi*X;
    %X2(:,1) = SP(Y(:,1), Phi, 8);
    r = zeros(nlabs,1);  %residual to each class
    %r2 = zeros(nlabs,1);
    for j = 1: nlabs
        index = find(D_lab(1,:) == labcls(j)); 
        A = D_lab(2:size(D_lab,1), index);
        %A = Phi(1:19,index);
        r(j) = norm(A*X(index,1)-Y(:,1),2);
       % r2(j) = norm(A*X2(index,1)-Y(:,1),2);
    end
    [~,tmp] = min(r);
    OMP_pred = labcls(tmp);
    %[~,SP_pred]=min(r2);
    
%end
%acc = length(find((facecls(test) - pred)==0))/length(test);






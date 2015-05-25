rho = [0.6 0.4]; B = [0.7 0.62; 0.62 0.75]; n1 = 400; 
graphMat1 = zeros(400, n1);
for i = 1: n1
    [A, ~, ~] = sbm(rho, B, 20);
    graphMat1(:, i) = A(:);
end

rho = [0.6 0.4]; B = [0.4 0.2; 0.2 0.4]; n2= 400;
graphMat2 = zeros(400, n2);
for i = 1: n2
    [A, ~, ~] = sbm(rho, B, 20);
    graphMat2(:, i) = A(:);
end

graphMat = [graphMat1 graphMat2];
graphLab = [ones(n1,1); ones(n2,1)+1];

graphInd = 1:(n1+n2);
ngraph = length(graphInd);

pred = zeros();
for i = 1:ngraph
    test = i;
    train = graphInd(setdiff(1:ngraph,test));
    test_data = graphMat(:, test);
    train_data = graphMat(:, train);
    [pred(i), ~] =  sparse_rep(test, train, test_data, train_data, graphLab); 
end

err = mean(pred' ~= graphLab);
function [ err, testResult ] = srcRepErrFun( Tau, Adj, s )


Tau= double(Tau);
[d,~]=size(Adj);
pred = zeros(d,1);
testResult = cell(d, 4);
tic
for i = 1:d
    
    tmp = Adj(1:d ~= i, 1:d ~= i);
    y = Adj(:,i);
    y = y(1:d ~= i);
    [testResult{i, 1}, testResult{i, 2}, testResult{i, 3}, testResult{i, 4}] = srcRep(i, y, tmp, Tau, s);
    pred(i) = testResult{i, 1};
end
toc;
%elapsedTime = toc;
%err = sum(abs(pred(:,1)-Tau))/d;
err = mean(pred ~= Tau);

end

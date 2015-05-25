function [ err, elapsedTime ] = errfun( Tau, Adj, s )


Tau= double(Tau);
[d,~]=size(Adj);
pred= zeros(d,1);
start=tic;
for i = 1:d
    
    tmp = Adj(1:d ~= i, 1:d ~= i);
    y = Adj(:,i);
    y = y(1:d ~= i);
    [pred(i),~] = graph_src(i, y, tmp, Tau, s);
    
end
elapsedTime = toc(start);
%err = sum(abs(pred(:,1)-Tau))/d;
err = mean(pred ~= Tau);

end

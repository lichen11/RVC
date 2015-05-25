function err = errfun( Tau, Adj, s )


Tau= double(Tau);
[d,~]=size(Adj);
pred= zeros(d,1);
%tic
for i = 1:d
    
    tmp = Adj(1:d ~= i, 1:d ~= i);
    y = Adj(:,i);
    y = y(1:d ~= i);
    [pred(i),~] = graph_src(i, y, tmp, Tau, s);
    
end
%toc;
%elapsedTime = toc;
%err = sum(abs(pred(:,1)-Tau))/d;
err = mean(pred ~= Tau);

end

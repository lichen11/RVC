function [err, elapsedTime ] = errsvdlda(tau, embed )
tau=double(tau);
nver=length(tau);
pred=zeros(nver,1);
start=tic;
for j=1:nver
    
    pred(j) = svdlda(j, embed, tau);
    
end
elapsedTime = toc(start);
err = mean(pred ~= tau);

end



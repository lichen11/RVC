function [ err, elapsedTime ] = errsvdknn( tau, embed, k)

tau=double(tau);
nver=length(tau);
pred=zeros(nver,1);
start=tic;
for j=1:nver
    
    pred(j) = svdknn(j, embed, tau, k);
    
end
elapsedTime = toc(start);
%err = sum(abs(pred(:,1)-tau))/nver;
err = mean(pred ~= tau);

end



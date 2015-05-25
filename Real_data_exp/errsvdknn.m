function err = errsvdknn( tau, embed, k)

tau=double(tau);
nver=length(tau);
pred=zeros(nver,1);
%tic;
for j=1:nver
    
    pred(j) = svdknn(j, embed, tau, k);
    
end
%toc;
%elapsedTime = toc;
%err = sum(abs(pred(:,1)-tau))/nver;
err = mean(pred ~= tau);

end



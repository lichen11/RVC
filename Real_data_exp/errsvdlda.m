function err = errsvdlda(tau, embed )
tau=double(tau);
nver=length(tau);
pred=zeros(nver,1);
%tic;
for j=1:nver
    
    pred(j) = svdlda(j, embed, tau);
    
end
%toc;
%elapsedTime = toc;
err = mean(pred ~= tau);

%end



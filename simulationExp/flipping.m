function [err2] = flipping(p, adj, tau, d, s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%rho = [0.6 0.4]; %B = [0.7 0.32; 0.32 0.75]; n=150;
%[adj, tau,~]=sbm(rho, B, n);
%embed = svdembed(2,adj);


%p = 0.1:0.05:0.95;
err2=zeros(length(p),3);
%elap2=zeros(length(p),3);
[dim1,~]=size(adj);

%Cor_gr = zeros(dim1, dim1, length(p));

for k= 1:length(p)
  inc = round(dim1*p(k));
  %randn('seed', 8);
  i = randsample(1:dim1, inc);

  patch = adj; 
  %patch(i,i) = 0;  %occlusion
  patch(i,i)= 1-patch(i,i); %errorful: flip
  %Cor_gr(:,:,k) = patch;
  %embed = svdembed(2,patch);
  embed = adjSpecEmbed(patch, d);
  [err2(k,1), ~] = errfun(tau, patch, s);
  [err2(k,2), ~] = errsvdknn(tau, embed, 1);  
  [err2(k,3), ~] = errsvdlda(tau, embed);
end

end


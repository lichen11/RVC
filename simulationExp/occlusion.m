% input p: occlusion vector
% adjacency matrix from the original graph
% tau the membership
% d the embedded dimension for ASE
% s the sparsity level for SRC
% K the number of nearest neighbors
function [err2] = occlusion(p, adj, tau, d, s, K)
%Input: Occlusion rate: p. e.g p = 0.1:0.05:0.95;

%Output: Occlusion misclassification error
%
%%
%
%
err2=zeros(length(p),3);
%elap2=zeros(length(p),3);
[dim1,~]=size(adj);

%Cor_gr = zeros(dim1, dim1, length(p));

for k= 1:length(p)
  inc = round(dim1*p(k));
  %seed = rng;
  i = randsample(1:dim1, inc);

  patch = adj; 
  patch(i,i) = 0;  %occlusion
  %patch(i,i)= 1-patch(i,i); %errorful: flip
  %Cor_gr(:,:,k) = patch;
  %embed = svdembed(2,patch);
  embed = adjSpecEmbed(patch, d);
  [err2(k,1), ~] = errfun(tau, patch, s);
  [err2(k,2), ~] = errsvdknn(tau, embed, K);  
  [err2(k,3), ~] = errsvdlda(tau, embed);
end

end

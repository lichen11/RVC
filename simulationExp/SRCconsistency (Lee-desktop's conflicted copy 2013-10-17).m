% sample the probability matrix from the B matrix
B = [ 0.7 0.35; 0.35 0.75];
rho = [0.6 0.4];
n = 2000;
[P, tau] = getProbMat(B,  rho, n);
A = rgSample(P);
ASE = adjSpecEmbed(A, 2);
predASE = kmeans(ASE, 2);
adjrand(predASE, tau)

[coeff1, residual1] = srcclust(P,1);
[coeff2, residual2] = srcclust(A*A,1);

pred = gmdistribution.fit(coeff1, 2);
%pred = kmeans(coeff1, 2);
adjrand(pred, tau)
pred2 = kmeans(coeff2, 2);
adjrand(pred2, tau) %  0.9349

class1 = find(tau==1);
class2 = find(tau==2);

mean(result(class1)) - mean(result(class2))
mean(result2(class1)) - mean(result2(class2))


tmp(tmp==0)=[]
A = reshape(A,4,3);


n = 500:100:1900
%for n_vec 





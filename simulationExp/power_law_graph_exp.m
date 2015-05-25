warning ('off','all');
N = 200;
M = 1;
alpha = 1.5%1.2;
b = 3%3;
X = randp(N,M,alpha,b);
X = ceil(X);
A = graph_from_degree_sequence(X);
tau = (X>5) + 1;

% embed = adjSpecEmbed(A, 10);
% err = zeros(10,3); 
% for d = 1:10
%     err(d,1) = errfun(tau, A, d);
%     err(d,2) = errsvdknn(tau, embed(:,1:d), 1);
%     err(d,3) = errsvdlda(tau, embed(:,1:d));
% end

d_vec = 1:140;
parfor i=1:length(d_vec)
    d = d_vec(i);
    s = d;
    embed = adjSpecEmbed(A, d);

    tmp1 = errfun(tau, A, s);
    tmp2 = errsvdknn(tau, embed, 1);
    tmp3 = errsvdlda(tau, embed);
    err(i,:) = [tmp1 tmp2 tmp3];
end








B = [0.7 0.32; 0.32 0.75]; 
n = 200;
rho = [0.4 0.6];
err2=zeros(1,3);

p_perm = randperm(n);
tau = [ones(1, n*rho(1)) ones(1, n*rho(2))+1];
off_diag_block = normrnd(0.32, 0.1, n*rho(1), n*rho(2));
diag_block1 = normrnd(0.7, 0.1, n*rho(1), n*rho(1));
diag_block2 = normrnd(0.75, 0.05, n*rho(2), n*rho(2));
P = [diag_block1 off_diag_block; off_diag_block' diag_block2];

tau = tau(p_perm)';
P = P(p_perm, p_perm);
U = rand(n);
U = triu(U,1);
U = U + U';
A = (U<P)+0;
diagIdx = 1:n+1:n^2;
A(diagIdx)=0;

embed = adjSpecEmbed(A, 10);
 
for d = 1:10
     errfun(tau, A, d)
     errsvdknn(tau, embed(:,1:d), 1)
     errsvdlda(tau, embed(:,1:d))
end
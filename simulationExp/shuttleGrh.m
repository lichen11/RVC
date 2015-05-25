%shuffle graph
%randomly shuffle indices
n = 25;  rho = [0.4 0.6]; B = [ 0.7 0.3; 0.3 0.75];

[A, tau, ~] = sbm(rho, B, n);
shuffle_ind = randperm(floor(n/2));
entire_ind = [shuffle_ind (floor(n/2)+1):n];
A_shuffle = A(entire_ind, entire_ind);
tau_shuffle = tau(entire_ind);

imagesc(A_shuffle)

[err1, testResult1] = srcRepErrFun(tau, dct2(A), 5);
[err2, testResult2] = srcRepErrFun(tau_shuffle, dct2(A_shuffle), 5);

mismatchRate1 = zeros(n, 1);




for i = 1 : n
   mismatchRate1(i) = 1 - sum(ismember(testResult1{entire_ind(i), 4}, testResult2{i, 4}))/length(testResult1{i,2});

end





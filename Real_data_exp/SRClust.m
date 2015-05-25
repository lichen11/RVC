%misLabRatio = zeros();
%for i = 1:50
rho = [0.6 0.4]; B = [0.7 0.32; 0.32 0.75];

[A, tau, ~] = sbm(rho, B, 150);
s = 75;
[err, result] = srcRepErrFun(tau, A, s);
X = reshape(cell2mat(result(:,4)), s, length(cell2mat(result(:,4)))/s);
X = X';
options = statset('Display','final');

gFitRaw = gmdistribution.fit(X, 2, 'options', options);

clsRaw = cluster(gFitRaw, X);
misLabRatioRaw = mean(clsRaw ~= tau);
adjrand = adjrand(clsRaw, tau)

D = pdist(X,'euclidean');

[Y,e] = cmdscale(D);

gFit = gmdistribution.fit(Y(:, 1:2), 2, 'options', options);
cls = cluster(gFit, Y(,1:2));
misLabRatio = mean(cls ~= tau);
adjrand = adjrand(tau, tau)
%end

embed = svdembed(2, A);
gFit2 = gmdistribution.fit(embed, 2, 'options', options);
cls2 = cluster(gFit2, embed);
misLabRatio2 = mean(cls2 ~= tau);


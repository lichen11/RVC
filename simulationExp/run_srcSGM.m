

nonZeroInd = find( sum(A) ~= 0);
A(nonZeroInd, nonZeroInd)

nonZeroIndAchem = find( sum(Achem) ~= 0);

[err, result] = srcRepErrFun(label(nonZeroInd), full(A(nonZeroInd, nonZeroInd)), 50);
[errAchem, ResultAchem] = srcRepErrFun(label(nonZeroIndAchem), full(binarA(Achem(nonZeroIndAchem, nonZeroIndAchem))), 50);

mismatchRate1 = zeros(length(nonZeroInd), 1);
mismatchRate2 = zeros(length(nonZeroInd), 1);
for i = 1 : length(nonZeroInd)
   mismatchRate1(i) = 1 - sum(ismember(result{i, 2}, ResultAchem{i, 2}))/length(result{i,2});
   mismatchRate2(i) = 1 - sum(ismember(ResultAchem{i, 2}, result{i ,2}))/length(ResultAchem{i, 2});
end

figure; plot(nonZeroInd, mismatchRate1); hold on; plot(nonZeroInd, mismatchRate2, 'r'); hold off;


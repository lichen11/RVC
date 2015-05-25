clc; clear all;
load herm_connectomes; load label; clear Achem Neuron_ordered;

[Adj, tau] = sortAdj(Agap, label);
Adj(logical(eye(size(Adj)))) = 0;
Adj=Adj~=0;
Adj = double(Adj);

Adj = full(Adj);

[err, elap] = errfun(label, Adj);
[ err, elap ] = errfun(label, full(Agap))
graph_src(5, Adj(:,5), Adj(1:279 ~= 5, 1:279 ~= 5), label)

embed = svdembed(3, full(Adj));
err = zeros(3,1);
[err(1), ~] = errfun(tau, Adj);
[err(2), ~] = errsvdknn(tau, embed, 1);
[err(3), ~] = errsvdlda(tau, embed);

embed2 = svdembed(3, full(Agap));
err2 = zeros(3, 1);
[ err2(1), elap ] = errfun(label, full(Agap));
[err2(2), ~] = errsvdknn(label, embed2, 1);
[err2(3), ~] = errsvdlda(label, embed2);





































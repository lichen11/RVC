tic;
p = 0.1:0.05:0.25;
err2=zeros(length(p),1);
elap2=zeros(length(p),1);
[dim1,~]=size(Adj);

Cor_gr = zeros(dim1, dim1, length(p));

for k= 1:length(p)
  inc = round(dim1*p(k));
  i = randsample(1:dim1, inc);

  patch = Adj; 
  patch(i,i) = 0;  %occlusion
  %patch(i,i)= 1-patch(i,i); %errorful: flip
  Cor_gr(:,:,k) = patch;
  embed = svdembed(2,patch);
  [err2(k), elap2(k)] = errfun(label, patch);
  %[err2(k,2), elap2(k,2)] = errsvdknn(tau, embed, 1);  
  %[err2(k,3), elap2(k,3)] = errsvdlda(tau, embed);
end

toc








rho = [0.6 0.4]; B = [0.7 0.32; 0.32 0.75]; %n=150;

err1 = occlusion(B, 145);
err2 = flipping(B, 145);
plot(0.1:0.05:0.95, err1,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest');
xlabel('Occlusion rate','fontsize',14)
ylabel('Error rate','fontsize',14)
figure;
plot(0.1:0.05:0.95, err2,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
%%
tic;
err_tensor1 = zeros(length(.1:.05:.95), 3, 30);
err_tensor2 = zeros(length(.1:.05:.95), 3, 30);
for nsim = 1:30
    err1 = occlusion(B, 145);
    err2 = flipping(B,145);
    err_tensor1(:,:, nsim) = err1;
    err_tensor2(:, :, nsim) = err2;
end
toc

for i = 1:30
    subplot(5,6,i);
    plot(0.1:0.05:0.95, err_tensor1(:,:,i),'-');
end

for i = 1:30
    subplot(5,6,i);
    plot(0.1:0.05:0.95, err_tensor2(:,:,i),'-');
end

x = mean(err_tensor1,3);
y = std(err_tensor1,1, 3); 
shadedErrorBar(0.1:0.05:0.95, x(:,1), y(:,1),'b',1);
hold on;
shadedErrorBar(0.1:0.05:0.95, x(:,2), y(:,2),'g',1);
shadedErrorBar(0.1:0.05:0.95, x(:,3), y(:,3),'r',1);
hold off;
xlabel('Occlusion rate','fontsize',14)
ylabel('Error rate','fontsize',14)


x = mean(err_tensor2,3);
y = std(err_tensor2,1, 3); 
shadedErrorBar(0.1:0.05:0.95, x(:,1), y(:,1),'b',1);
hold on;
shadedErrorBar(0.1:0.05:0.95, x(:,2), y(:,2),'g',1);
shadedErrorBar(0.1:0.05:0.95, x(:,3), y(:,3),'r',1);
hold off;
xlabel('Flipping rate','fontsize',14)
ylabel('Error rate','fontsize',14)


%  shadedErrorBar(0.1:0.05:0.95, err_tensor1(:,1,1) ,{@mean,@std},'-r',1); 
hold on;
plot(0.1:0.05:0.95, mean(err_tensor1,3) + std(err_tensor1,1, 3), '--','LineWidth',...
    2,'MarkerSize',8 )
plot(0.1:0.05:0.95, mean(err_tensor1,3) - std(err_tensor1,1,3), '--','LineWidth',...
    2,'MarkerSize',8 )
plot([0, 1], [0.4, 0.4], 'y-')

plot(0.1:0.05:0.95, mean(err_tensor2,3), 'o-','LineWidth',2,'MarkerSize',8 )
hold on;
plot(0.1:0.05:0.95, mean(err_tensor2,3) + std(err_tensor2,1, 3), '--','LineWidth',...
    2,'MarkerSize',8 )
plot(0.1:0.05:0.95, mean(err_tensor2,3) - std(err_tensor2,1,3), '--','LineWidth',...
    2,'MarkerSize',8 )
plot([0, 1], [0.4, 0.4], 'y-')
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
%xlabel('Occlusion rate','fontsize',14)
xlabel('Flipping rate','fontsize',14)
ylabel('Error rate','fontsize',14)
%x = mean(err_tensor1, 3);
x = mean(err_tensor2, 3);
tmp = 0.1:0.05:0.95;
plot(tmp, x(:,1), '-*k',tmp, x(:,2), '-ok',tmp, x(:,3),...
    '-^k', 'LineWidth', 2, 'MarkerSize',8)

%plot(0.1:0.05:0.95, mean(err_tensor1,3), 'o-','LineWidth',2,'MarkerSize',8 )
%plot(0.1:0.05:0.95, mean(err_tensor2,3), 'o-','LineWidth',2,'MarkerSize',8 )
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest')
xlabel('Occlusion rate','fontsize',14)
xlabel('Reversion rate','fontsize',14)
ylabel('Error rate','fontsize',14)
savefile = 'err_tensor.mat';
save(savefile, 'err_tensor1', 'err_tensor2')
%%
n = 15:10:250;
err = zeros(length(n),3);
for i = 1:length(n)
   [adj, tau, ~] = sbm(rho, B, n(i));
   embed = svdembed(2, adj);
   [err(i,1), ~] = errfun(tau, adj);
   [err(i,2), ~] = errsvdknn(tau, embed, 1);  
   [err(i,3), ~] = errsvdlda(tau, embed);
end

plot(n, err(:,1),'-*k', n, err(:,2),  '-ok', ...
    n, err(:, 3),  '-^k','LineWidth',2,'MarkerSize',8);
hold on;
line(n,[0.4 0.4]); 
legend('SRC','NN o SVD','LDA o SVD','Location','SouthWest');
xlabel('Vertices','fontsize',14)
ylabel('Error rate','fontsize',14)
[adj, tau,~]=sbm(rho, B,n);
%embed = svdembed(2,adj);


p = 0.1:0.05:0.95;
err2=zeros(length(p),3);
elap2=zeros(length(p),3);
[dim1,~]=size(adj);

Cor_gr = zeros(dim1, dim1, length(p));

for k= 1:length(p)
  inc = round(dim1*p(k));
  i = randsample(1:dim1, inc);

  patch = adj; 
  %patch(i,i) = 0;  %occlusion
  patch(i,i)= 1-patch(i,i); %errorful: flip
  Cor_gr(:,:,k) = patch;
  embed = svdembed(2,patch);
  [err2(k,1), elap2(k,1)] = errfun(tau, patch);
  [err2(k,2), elap2(k,2)] = errsvdknn(tau, embed, 1);  
  [err2(k,3), elap2(k,3)] = errsvdlda(tau, embed);
end
figure;
plot(0.1:0.05:0.95, err2,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
xlabel('Flipping rate','fontsize',14)
ylabel('Error rate','fontsize',14)

%%





%Clean graph data
incr = 15:10:145;
err=zeros(length(incr),3);
elap = zeros(length(incr),3);
for i=1:length(incr)
   stemp = sprintf('SBM%d',incr(i));
   load(stemp);
   [err(i,1), elap(i,1)] = errfun(tau, adj);
   [err(i,2), elap(i,2)] = errsvdknn(tau, adj, 2, 1);  
   [err(i,3), elap(i,3)] = errsvdlda(tau, adj, 2);
      
end


plot(incr, err,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','SVD o KNN','SVD o LDA','Location','SouthWest')
xlabel('Vertices','fontsize',14)
ylabel('Error rate','fontsize',14)
figure;
plot(incr, elap,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','SVD o KNN','SVD o LDA','Location','SouthWest')
xlabel('Vertices','fontsize',14)
ylabel('Computation time','fontsize',14)
%%
%corrupted data

p = 0.1:0.05:0.85;
err2=zeros(length(p),3);
elap2=zeros(length(p),3);
[dim1,~]=size(adj);

Cor_gr = zeros(dim1, dim1, length(p));

for k= 1:length(p)
  inc = round(dim1*p(k));
  i = randsample(1:dim1, inc);

  patch = adj; 
  %patch(i,i) = 0;  %occlusion
  patch(i,i)= 1-patch(i,i); %errorful: flip
  Cor_gr(:,:,k) = patch;
  
  [err2(k,1), elap2(k,1)] = errfun(tau, patch);
  [err2(k,2), elap2(k,2)] = errsvdknn(tau, patch, 2, 1);  
  [err2(k,3), elap2(k,3)] = errsvdlda(tau, patch, 2);
end

figure;
plot(p, err2,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','SVD o KNN','SVD o LDA','Location','SouthWest')
xlabel('Contamination percentage','fontsize',14)
ylabel('Error rate','fontsize',14)
figure;
plot(p, elap2,'LineWidth',2,'MarkerSize',8);
legend('SRC','SVD o KNN','SVD o LDA','Location','SouthWest')
xlabel('Contamination percentage','fontsize',14)
ylabel('Computation time','fontsize',14)
%%
%averaging
p = 0.1:0.05:0.85;
err_tensor = zeros(length(p), 3, 10);
elap_tensor = zeros(length(p), 3, 10);
for nsim=1:10
   stemp = sprintf('n145SBM%d',nsim);
   load(stemp);

    
    err3=zeros(length(p),3);
    elap3=zeros(length(p),3);
    [dim1,~]=size(adj);

    Cor_gr = zeros(dim1, dim1, length(p));

    for k= 1:length(p)
        inc = round(dim1*p(k));
        i = randsample(1:dim1, inc);

        patch = adj; 
        patch(i,i) = 0;
        Cor_gr(:,:,k) = patch;
  
        [err3(k,1), elap3(k,1)] = errfun(tau, patch);
        [err3(k,2), elap3(k,2)] = errsvdknn(tau, patch, 2, 1);  
        [err3(k,3), elap3(k,3)] = errsvdlda(tau, patch, 2);
    end

    err_tensor(:,:,nsim)=err3;
    elap_tensor(:,:,nsim)=elap3;


end
savefile='corr_result.mat'
save(savefile, 'err_tensor','elap_tensor')

for j=1:10
    subplot(2,5,j)
    plot(p, err_tensor(:,:,j))
   ylabel('Error rate')
   xlabel('missing percentage')
end

 legend('SRC','SVDoKNN','SVDoLDA','Location','SouthWest')
    
    

for j=1:9
    subplot(3,3,j)
    plot(p, elap_tensor(:,:,j))
    
end

 plot(p, mean(err_tensor,3), 'o-','LineWidth',2,'MarkerSize',8 )
 plot(p, mean(elap_tensor,3), 'o-','LineWidth',2,'MarkerSize',8 )
legend('SRC','SVD o KNN','SVD o LDA','Location','SouthWest')
    xlabel('Contamination percentage','fontsize',14)
    ylabel('Computing time','fontsize',14)
     ylabel('Error rate','fontsize',14)

 for j = 1:16    
    subplot(4,4,j)
    imshow(Cor_gr(:,:,j));
 end
%% 
incr = 0:0.05:0.48;
err=zeros(length(incr),3);
elap = zeros(length(incr),3);
for i=1:length(incr)
    
   stemp = sprintf('../homogenous_data/145one_paramSBM%f.mat',incr(i));
   load(stemp);
   [err(i,1), elap(i,1)] = errfun(tau, adj);
   [err(i,2), elap(i,2)] = errsvdknn(tau, adj, 2, 1);  
   [err(i,3), elap(i,3)] = errsvdlda(tau, adj, 2);
      
end
figure;
plot(incr, err,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
xlabel('Number of vertices','fontsize',14)
ylabel('Error rate','fontsize',14)
%set(gca,'XTickLabel',incr); 


 
 
%% 
incr = 15:10:145;
err=zeros(length(incr),3);
elap = zeros(length(incr),3);
for i=1:length(incr)
   stemp = sprintf('../homogenous_data/denseSBM%d.mat',incr(i));
   load(stemp);
   
   
   [err(i,1), elap(i,1)] = errfun(tau, adj);
   [err(i,2), elap(i,2)] = errsvdknn(tau, adj, 2, 1);  
   [err(i,3), elap(i,3)] = errsvdlda(tau, adj, 2);
      
end
figure;
plot(incr, err,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
xlabel('alpha','fontsize',14)
ylabel('Error rate','fontsize',14)

%%
p = 0.1:0.05:0.85;
incr = 15:10:145;
err3dmat = zeros(length(p),4,length(incr));

err=zeros(length(incr),3);
elap = zeros(length(incr),3);
for i=1:length(incr)
   stemp = sprintf('../homogenous_data/denseSBM%d.mat',incr(i));
   load(stemp);
  
err2=zeros(length(p),3);
elap2=zeros(length(p),3);
[dim1,~]=size(adj);

Cor_gr = zeros(dim1, dim1, length(p));

for k= 1:length(p)
  inc = round(dim1*p(k));
  j = randsample(1:dim1, inc);

  patch = adj; 
  %patch(i,i) = 0;  %occlusion
  patch(j,j)= 1-patch(j,j); %errorful: flip
  Cor_gr(:,:,k) = patch;
  
  [err2(k,1), elap2(k,1)] = errfun(tau, patch);
  [err2(k,2), elap2(k,2)] = errsvdknn(tau, patch, 2, 1);  
  [err2(k,3), elap2(k,3)] = errsvdlda(tau, patch, 2);
end
   
   %figure;
   err3dmat(:,:,i) = [err2 p'];
   %[err(i,1), elap(i,1)] = errfun(tau, adj);
   %[err(i,2), elap(i,2)] = errsvdknn(tau, adj, 2, 1);  
   %[err(i,3), elap(i,3)] = errsvdlda(tau, adj, 2);
      
end
for j = 1:14
subplot(3,5,j)
plot(p, err3dmat(:,1:3,j))
end 
plot(p, incr, err3dmat(:,1:3,:))
%%
rho = [0.4 0.6];
B = [0.7 0.3; 0.3 0.75];
incr = 15:10:145;
err=zeros(length(incr),3);
elap = zeros(length(incr),3);
for i=1:length(incr)
   [adj, tau, ~] = sbm(rho,B,incr(i));
   [err(i,1), elap(i,1)] = errfun(tau, adj);
   [err(i,2), elap(i,2)] = errsvdknn(tau, adj, 2, 1);  
   [err(i,3), elap(i,3)] = errsvdlda(tau, adj, 2);
      
end
figure;
plot(incr, err,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
xlabel('Number of vertices','fontsize',14)
ylabel('Error rate','fontsize',14)
%%
rho = [0.4 0.6];
incr = 0:0.05:0.45;
err=zeros(length(incr),3);
elap = zeros(length(incr),3);
for i=1:length(incr)
   B = [0.4+incr(i) 0.3; 0.3 0.45+incr(i)];
   [adj, tau, ~] = sbm(rho,B,200);
   
   [err(i,1), elap(i,1)] = errfun(tau, adj);
   [err(i,2), elap(i,2)] = errsvdknn(tau, adj, 2, 1);  
   [err(i,3), elap(i,3)] = errsvdlda(tau, adj, 2);
      
end
figure;
plot(incr, err,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
xlabel('alpha','fontsize',14)
ylabel('Error rate','fontsize',14)
%%
%randn('seed',0)
p = 0.1:0.05:0.95;
rho = [0.4 0.6];
incr = 0:0.05:0.45;
%elap = zeros(length(incr),3);
err3dmat = zeros(length(p),4,length(incr));
for i=1:length(incr)
   B = [0.4+incr(i) 0.3; 0.3 0.45+incr(i)];
   [adj, tau, ~] = sbm(rho,B,150);
   [dim1,~]=size(adj);
   err=zeros(length(incr),3);
   Cor_gr = zeros(dim1, dim1, length(p));
   for k= 1:length(p)
        inc = round(dim1*p(k));
        j = randsample(1:dim1, inc);
        patch = adj; 
        patch(j,j) = 0;  %occlusion
        %patch(i,i)= 1-patch(i,i); %errorful: flip
        Cor_gr(:,:,k) = patch;
 
        [err(k,1), ~] = errfun(tau, adj);
        [err(k,2), ~] = errsvdknn(tau, adj, 2, 1);  
        [err(k,3), ~] = errsvdlda(tau, adj, 2);
   end
   err3dmat(:,:,i) = [err p'];   
end
figure;
plot(incr, err,'o-','LineWidth',2,'MarkerSize',8);
legend('SRC','KNN o SVD','LDA o SVD','Location','SouthWest')
xlabel('alpha','fontsize',14)
ylabel('Error rate','fontsize',14)



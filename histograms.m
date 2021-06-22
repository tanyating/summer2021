addpath('utils');
clf;

Rs = [0, pi/2, pi, pi*3/2]; %4 rotations
N = 100; %N grids
p = 10; %molecule length (p>q)
q = 10; %molecule width
mol = rand(q,p).*randi(10,q,p); %random molecules in 2D
ts = 0:N-p; %N-p+1 translations
A = zeros(N-p+1,4,N); %store all a_{t,R} (index:translation, rotation, vector)

%construct a_{t,R} based on mol
for i=0:N-p
    for j=1:4
        tmp=sum(rot90(mol,j-1),1);
        A(1+i,j,1+i:i+length(tmp)) = tmp;
    end
end

M = 10000; %number of noise vectors
sigma = 10; %level of noise
cov = sigma^2.*eye(N);
x = mvnrnd(zeros(N,1), cov, M); %iid noise
xhat = (x-mean(x,2))./sqrt(sum((x-mean(x,2)).^2,2)); %xhat
d1s = zeros(M,N-p+1,4); %metric d1: |x-a_{t,R}|
d2s = zeros(M,N-p+1,4); %metric d1: <xhat,a_{t,R}hat>
fs = zeros(M,N-p+1); %method 3: f(t) = max_R <xhat,a_{t,R}hat>
js = zeros(M,N-p+1); %method 4

a = zeros(1,N);
ahat = zeros(1,N);
for i=0:N-p
   cur_prod = zeros(M,4);
   for j=1:4
       a(1,:) = A(i+1,j,:);
       ahat = (a-mean(a))/norm(a-mean(a));
       d1s(:,i+1,j) = sum((x-a).^2,2); %distance to current a_{t,R}
       d2s(:,i+1,j) = xhat*transpose(ahat); %inner product <xhat,a_{t,R}hat>
       cur_prod(:,j) = xhat*transpose(ahat);
   end
   fs(:,i+1) = max(cur_prod,[],2); %f(t=i)
   js(:,i+1) = (fs(:,i+1)-mean(cur_prod,2))./std(cur_prod,0,2); %j(t=i)
end

nbins=100;
k=1;
figure(k);
histogram(d1s,nbins);
title('d1');
k=k+1;

figure(k);
for i=1:4
    histogram(d1s(:,:,i),nbins);
    hold on;
end
legend('R1', 'R2', 'R3', 'R4');
title('d1');
k=k+1;

figure(k);
% binwidth = (max(d2s,[],'all')-min(d2s,[],'all'))/nbins;
histogram(d2s,nbins,'Normalization','probability');
title('d2');
k=k+1;

figure(k);
for i=1:4
    histogram(d2s(:,:,i),nbins);
    hold on;
end
legend('R1', 'R2', 'R3', 'R4');
title('d2');
k=k+1;

figure(k);
histogram(fs,nbins);
title('f(t)');
k=k+1;

figure(k);
histogram(js,nbins);
title('j(t)');


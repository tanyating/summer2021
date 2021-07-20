clear;
addpath('utils','dp');


N = 24; %N grids = # pixels in 1D
p = 4; %molecule length (p>=q)
q = 4; %1 %molecule width
Nc = (N-p+1)*4; % number of configurations (signal)
Nt = N-p+1; % number of translations
mol = molecule(p,q) %.*randi(10,q,p); %random molecule in 2D

%construct a_{t,R} based on mol
A = template(mol,N); 

M = 20000; %number of noise vectors
sigma = 0.1; %level of noise
cov = sigma^2.*eye(N);
p_0 = 0.5; % prior prob for noise (no signal)
[y,tl_class] = randdata(M,A,sigma,p_0); % generate y and true labels
x = y(tl_class==0,:); % iid noise
y = y(tl_class~=0,:); % positive (signal) instances

M1 = length(x); % number of noise
M2 = length(y); % number of positive instances

x_d1s = zeros(M1,N-p+1,4); % d1: |x-a_{t,R}| for noise x 
x_d2s = zeros(M1,N-p+1,4); % d2: <xhat,a_{t,R}hat> for noise x
x_fs = zeros(M1,N-p+1); % method 3: f(t) = max_R <xhat,a_{t,R}hat> for noise x
x_js = zeros(M1,N-p+1); % method 4 for noise x

y_fs = zeros(M2,N-p+1); % method 3: f(t) for signals y

a = zeros(1,N);
ahat = zeros(1,N);
for i=0:N-p
   x_cur_prod = zeros(M1,4);
   y_cur_prod = zeros(M2,4);
   for j=1:4
       k=map_class([i j]);
       a(1,:) = A(k,:);
       x_d1s(:,i+1,j) = d1(x,a); %distance to current a_{t,R}
       x_d2s(:,i+1,j) = d2(x,a); %inner product <xhat,a_{t,R}hat>
       x_cur_prod(:,j) = d2(x,a);
       y_cur_prod(:,j) = d2(y,a);
   end
   x_fs(:,i+1) = max(x_cur_prod,[],2); %f(t=i)
   y_fs(:,i+1) = max(y_cur_prod,[],2);
   x_js(:,i+1) = (x_fs(:,i+1)-mean(x_cur_prod,2))./std(x_cur_prod,0,2); %j(t=i)
end

nbins=100;

% figure;
% histogram(x_d1s,nbins);
% title('d1');
% 
% figure;
% for i=1:4
%     histogram(x_d1s(:,:,i),nbins);
%     hold on;
% end
% legend('R1', 'R2', 'R3', 'R4');
% title('d1');
% 
% figure;
% % binwidth = (max(d2s,[],'all')-min(d2s,[],'all'))/nbins;
% histogram(x_d2s,nbins,'Normalization','probability');
% title('d2');
% 
% figure;
% for i=1:4
%     histogram(x_d2s(:,:,i),nbins);
%     hold on;
% end
% legend('R1', 'R2', 'R3', 'R4');
% title('d2');
% 
% figure;
% histogram(x_fs,nbins);
% title('f(t)');
% 
% figure;
% histogram(x_js,nbins);
% title('j(t)');

figure;
histogram(max(x_fs,[],2),nbins);
hold on;
histogram(max(y_fs,[],2),nbins);
title('max inner product');
legend('no signal', 'signal');



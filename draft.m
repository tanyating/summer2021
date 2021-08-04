clear;
addpath('utils','dp');

p=8;
q=32; % length (number of nonzero entries on N pixels)
sigma=5;
f = @(xx)(exp(-xx.^2/(2.*p)).*erf(xx/(2*sqrt(2)*sigma)).*xx.^(q-1));


A(1) = 0;
A(2) = 2;
A(3) = 2*pi;
for i=4:q+1
    A(i) = 2*pi/(i-2)*A(i-2);
end

efpr = 1/2 - 1/2*A(q+1)*1/(sqrt(2*pi*p).^q)*integral(f,0,inf)

%%% check for equal dist y1 for different R
%         ds_cur = ds_cur.^2;
%         d1_cur = ds_cur(:,1:4:end);
%         d1s(j,i,:) = d1_cur(:);
%         d2_cur = ds_cur(:,2:4:end);
%         d2s(j,i,:) = d2_cur(:);
%         d3_cur = ds_cur(:,3:4:end);
%         d3s(j,i,:) = d3_cur(:);
%         d4_cur = ds_cur(:,4:4:end);
%         d4s(j,i,:) = d4_cur(:);
% mean(mean(d1s,3),2)
% mean(mean(d2s,3),2)

%%% naming for stacked bars
% xticks(1:length(ps));
% names = strings(length(ps));
% for i=1:length(ps)
%     names(i) = sprintf('p=%d q=%d', ps(i),qs(i));
% end
% xticklabels(names);


%%% original histograms
% clear;
% addpath('utils','dp');
% 
% raw = 0;   % AHB added: 0: use best dist as in meth2,   1: use all dists (Manas) giving tail as in Rickgauer Fig 1e.
% 
% seed =  0; rng(seed);
% 
% N = 24; %N grids = # pixels in 1D
% p = 8; %molecule length (p>=q)
% q = 2; %1 %molecule width
% Nc = (N-p+1)*4; % number of configurations (signal)
% Nt = N-p+1; % number of translations
% mol = molecule(p,q) %.*randi(10,q,p); %random molecule in 2D
% 
% %construct a_{t,R} based on mol
% A = template(mol,N); 
% 
% M = 20000; %number of noise vectors
% sigma = 1; %level of noise
% cov = sigma^2.*eye(N);
% p_0 = 0.5; % prior prob for noise (no signal)
% [y,tl_class] = randdata(M,A,sigma,p_0); % generate y and true labels
% x = y(tl_class==0,:); % iid noise
% y = y(tl_class~=0,:); % positive (signal) instances
% 
% M1 = length(x); % number of noise
% M2 = length(y); % number of positive instances
% 
% x_d1s = zeros(M1,N-p+1,4); % d1: |x-a_{t,R}| for noise x 
% x_d2s = zeros(M1,N-p+1,4); % d2: <xhat,a_{t,R}hat> for noise x
% 
% y_d2s = zeros(M2,N-p+1,4); % d2: <xhat,a_{t,R}hat> for noise x   for MANAS
% 
% x_fs = zeros(M1,N-p+1); % method 3: f(t) = max_R <xhat,a_{t,R}hat> for noise x
% x_js = zeros(M1,N-p+1); % method 4 for noise x
% 
% y_fs = zeros(M2,N-p+1); % method 3: f(t) for signals y
% 
% 
% a = zeros(1,N);
% ahat = zeros(1,N);
% for i=0:N-p
%    x_cur_prod = zeros(M1,4);
%    y_cur_prod = zeros(M2,4);
%    for j=1:4
%        k=map_class([i j]);
%        a(1,:) = A(k,:);
%        x_d1s(:,i+1,j) = d1(x,a); %distance to current a_{t,R}
%        x_d2s(:,i+1,j) = d2(x,a); %inner product <xhat,a_{t,R}hat>
%        y_d2s(:,i+1,j) = d2(y,a); %inner product <xhat,a_{t,R}hat>  for MANAS
%        x_cur_prod(:,j) = d2(x,a);
%        y_cur_prod(:,j) = d2(y,a);
%    end
%    x_fs(:,i+1) = max(x_cur_prod,[],2); %f(t=i)
%    y_fs(:,i+1) = max(y_cur_prod,[],2);
%    x_js(:,i+1) = (x_fs(:,i+1)-mean(x_cur_prod,2))./std(x_cur_prod,0,2); %j(t=i)
% end
% 
% nbins=100;

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

% figure;
% if raw
%   histogram(x_d2s(:),nbins); hold on;
%   histogram(y_d2s(:),nbins);
%   title('all inner products');
% else
%   histogram(max(x_fs,[],2),nbins); hold on;
%   histogram(max(y_fs,[],2),nbins);
%   title('max inner product');
% end
% legend('no signal', 'signal');


%%% avg over seed

%         fps = zeros(1,length(seeds));
%         tp1s = zeros(1,length(seeds));
%         tp2s = zeros(1,length(seeds));
%         tp3s = zeros(1,length(seeds));
%         
%         fat_tp1s = zeros(1,length(seeds));
%         tall_tp1s = zeros(1,length(seeds));
%         
%         for l=1:length(seeds)
%             
%             seed = seeds(l); % rand seed generator
%             rng(seed);
            
%         end
%         
%         fpmean(j,i) = mean(fps);
%         tp1mean(j,i) = mean(tp1s);
%         tp2mean(j,i) = mean(tp2s);
%         tp3mean(j,i) = mean(tp3s);

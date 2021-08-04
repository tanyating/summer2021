clear
addpath('utils','dp');

sigma = 5;
seed = 0;

N = 96; %N grids = # pixels in 1D
p = 32;%[10,8,7,4,9];
q = 8;%[2,2,3,4,3];
Nc = (N-p+1)*4; % number of configurations (signal)
Nt = N-p+1; % number of translations

tao = 0.31; % threshold for method 2

k=1;
for j=1:2 % iterate thru different methods
    

%     ds = zeros(1,length(seeds));
%     bs = zeros(1,length(seeds));
%     fps = zeros(1,length(seeds));
%     gs = zeros(1,length(seeds));
%     h1s = zeros(1,length(seeds));
%     h2s = zeros(1,length(seeds));
%     h3s = zeros(1,length(seeds));
%     h4s = zeros(1,length(seeds));
%     rs = zeros(1,length(seeds));
%     
%     for i=1:length(seeds) % iterate thru different seeds
%         
%         seed = seeds(i); % rand seed generator
        rng(seed);
        
        mol = molecule(p,q,seed); %random molecule in 2D
%         mol = randn(p,q);

        %construct a_{t,R} based on mol
        A = template(mol,N); 


        cov = sigma^2.*eye(N);
        
        M = 2000; % number of random examples
        p_0 = 0.5; % prior prob for noise (no signal)
        [y,tl_class] = randdata(M,A,sigma,p_0); % generate y and true labels
        tl_pairs = inverse_map(tl_class);
        
        % predict labels by minimizing distance (norm)
        if (j==1)
            [pl_class,ds_cur] = detect_max(y,A,@(y,a)-d1(y,a));
        else
            [pl_class,ds_cur] = detect_max(y,A,@(y,a)d2(y,a), tao);
        end
        
        pl_pairs = inverse_map(pl_class);
        
        fps(j) = sum(tl_class==0 & pl_class>0)/sum(tl_class==0); % false positive rate
        fat_fps(j) = sum(tl_class==0 & (pl_pairs(:,2)==1 | pl_pairs(:,2)==3))/sum(tl_class==0); % fat fp rate
        tall_fps(j) = sum(tl_class==0 & (pl_pairs(:,2)==2 | pl_pairs(:,2)==4))/sum(tl_class==0); % tall fp rate
        %             fns(k, l) = sum(tl_class>0 & pl_class==0)/sum(tl_class>0); % false negative rate
        
        
        tp1s(j) = sum(tl_class>0 & pl_class>0)/sum(tl_class>0); % overall tp rate
        fat_tp1s(j) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & pl_class>0)/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat tp rate
        tall_tp1s(j) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & pl_class>0)/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall tp rate
        
        tp2s(j) = sum(tl_class>0 & (pl_class==tl_class))/sum(tl_class>0); % correct (t,R) rate
        fat_tp2s(j) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & (pl_class==tl_class))/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat correct (t,R)
        tall_tp2s(j) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & (pl_class==tl_class))/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall correct (t,R)
        
        tp3s(j) = sum(tl_class>0 & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_class>0); % correct t rate
        fat_tp3s(j) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat correct (t,R)
        tall_tp3s(j) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall correct (t,R)
        
        %             image_show = abs(sigma-sigma_error)<1e-14;
        image_show = 0;
        
        %         C = error_matrix(tl_class,pl_class,Nc,image_show); % error matrix for (t,R) pair
        C_red = error_matrix_red((tl_class),(pl_class),Nc,Nt,image_show); % reduced error matrix for (t,R) pair
        
        [g,h1,h2,h3,h4,o1,o2,o3,o4,r] = extract_C(C_red,p,Nt,image_show);
        gs(j) = g; % avg rate for true t, wrong R per t
        h1s(j) = h1; % avg fp when R=1
        h2s(j) = h2; % avg fp when R=2
        h3s(j) = h3; % avg fp when R=3
        h4s(j) = h4; % avg fp when R=4
        rs(j) = r; % avg correct (t,R)
%         o1s(j, i) = o1; % approx fn when R=1
%         o2s(j, i) = o2; % approx fn when R=2
%         o3s(j, i) = o3; % approx fn when R=3
%         o4s(j, i) = o4; % approx fn when R=4
        
        %         Ct = error_matrix(get_tr(tl_class),get_tr(pl_class),Nt,image_show); % error matrix for translation t
        Ct_red = error_matrix_red(get_tr(tl_class),get_tr(pl_class),Nt,Nt,image_show); % error matrix for t
        
        % extract error rates from Ct_red
        [a,b,c,d,e,F] = extract_Ct(Ct_red,p);
        
%         cs(j, i) = c; % tn for noise
        ds(j) = d; % tp per translation
%         as(j, i) = a; % avg fp per translation
        bs(j) = b; % fn per translation
%         es(j, i) = e; % avg misclassifing rate (overlapping) per translation
%         Fs(j, i) = F; % noise rate (non-overlapping) per translation

%     end
%     
%     dmean(j) = mean(ds);
%     bmean(j) = mean(bs);
%     fpmean(j) = mean(fps);
%     h1mean(j) = mean(h1s);
%     h2mean(j) = mean(h2s);
%     h3mean(j) = mean(h3s);
%     h4mean(j) = mean(h4s);
%     rmean(j) = mean(rs);
%     gmean(j) = mean(gs);

end



figure;
hold on;
mol_bar = zeros(2,4);
mol_bar(:,1) = rs;
mol_bar(:,2) = ds-rs;
mol_bar(:,3) = 1-ds-bs;
mol_bar(:,4) = bs;
bar(1:2,mol_bar,'stacked');
names = strings(2);
for i=1:2
    names(i) = sprintf('method %d', i);
end
xticks(1:2);
xticklabels(names);
ylim([0 inf]);
title('actual mol');
legend('correct (t,R)', 'correct t, wrong R','wrong t', 'fn');

figure;
hold on;
noise_bar = zeros(2,2);
noise_bar(:,1) = 1-fps;
noise_bar(:,2) = fps;
bar(1:2,noise_bar,'stacked');
xticks(1:2);
xticklabels(names);
title('actual noise');
legend('tn','fp');

figure;
hold on;
rot_bar = zeros(2,2);
rot_bar(:,1) = h1s+h3s;
rot_bar(:,2) = h2s+h4s;
rot_bar = rot_bar./sum(rot_bar,2); % normalize over the 4 fps (sum to 1)
bar(1:2,rot_bar,'stacked');
xticks(1:2);
xticklabels(names);
legend('fat', 'tall');
title('False Positive rate per orientation');

figure;
hold on;
ratio_bar = zeros(2,1);
ratio_bar(:,1) = (h1s+h3s)./(h2s+h4s);
bar(1:2,ratio_bar,'stacked');
xticks(1:2);
xticklabels(names);
hline(1);
title('Ratio of fat FPR and tall FPR');




figure;
k=1;
% for j=1:2
ax = subplot(1,3,k);
hold on;

j=1;
plot(fps(j), tp1s(j),'b.', 'Markersize', 10);
plot(fat_fps(j), fat_tp1s(j),'r.', 'Markersize', 10);
plot(tall_fps(j), tall_tp1s(j),'g.', 'Markersize', 10);
j=2;
plot(fps(j), tp1s(j),'b+', 'Markersize', 10);
plot(fat_fps(j), fat_tp1s(j),'r+', 'Markersize', 10);
plot(tall_fps(j), tall_tp1s(j),'g+', 'Markersize', 10);

xx=0:0.01:1;
plot(xx,xx);
xlabel('general fp');
ylabel('general tp');
% xlabel(sprintf('p=%d q=%d sigma=%.2f', p,q,sigma));
% title('ROC (general fp, general tp)');

k=k+1;

ax = subplot(1,3,k);
hold on;

j=1;
plot(fps(j), tp2s(j),'b.', 'Markersize', 10);
plot(fat_fps(j), fat_tp2s(j),'r.', 'Markersize', 10);
plot(tall_fps(j), tall_tp2s(j),'g.', 'Markersize', 10);
j=2;
plot(fps(j), tp2s(j),'b+', 'Markersize', 10);
plot(fat_fps(j), fat_tp2s(j),'r+', 'Markersize', 10);
plot(tall_fps(j), tall_tp2s(j),'g+', 'Markersize', 10);

xx=0:0.01:1;
plot(xx,xx);
xlabel('general fp');
ylabel('correct (t,R)');
% xlabel(sprintf('p=%d q=%d sigma=%.2f', p,q,sigma));
% title('ROC (general fp, correct (t,R))');
k=k+1;

ax = subplot(1,3,k);
hold on;

j=1;
plot(fps(j), tp3s(j),'b.', 'Markersize', 10);
plot(fat_fps(j), fat_tp3s(j),'r.', 'Markersize', 10);
plot(tall_fps(j), tall_tp3s(j),'g.', 'Markersize', 10);
j=2;
plot(fps(j), tp3s(j),'b+', 'Markersize', 10);
plot(fat_fps(j), fat_tp3s(j),'r+', 'Markersize', 10);
plot(tall_fps(j), tall_tp3s(j),'g+', 'Markersize', 10);

xx=0:0.01:1;
plot(xx,xx);
xlabel('general fp');
ylabel('correct t');
% xlabel(sprintf('p=%d q=%d sigma=%.2f', p,q,sigma));
% title('ROC (general fp, correct t)');
k=k+1;
% end
legend('overall (method 1)','fat (method 1)','tall (method 1)','overall (method 2)','fat (method 2)','tall (method 2)');


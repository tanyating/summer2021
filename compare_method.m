clear
addpath('utils','dp');

sigma = 2;
seeds = 0:0;

N = 48; %N grids = # pixels in 1D
p = 16;%[10,8,7,4,9];
q = 4;%[2,2,3,4,3];
Nc = (N-p+1)*4; % number of configurations (signal)
Nt = N-p+1; % number of translations

tao = 0.44; % threshold for method 2

k=1;
for j=1:2 % iterate thru different methods
    

    ds = zeros(1,length(seeds));
    bs = zeros(1,length(seeds));
    fps = zeros(1,length(seeds));
    gs = zeros(1,length(seeds));
    h1s = zeros(1,length(seeds));
    h2s = zeros(1,length(seeds));
    h3s = zeros(1,length(seeds));
    h4s = zeros(1,length(seeds));
    rs = zeros(1,length(seeds));
    
    for i=1:length(seeds) % iterate thru different seeds
        
        seed = seeds(i); % rand seed generator
        rng(seed);
        
        mol = molecule(p,q,seed); %random molecule in 2D
%         mol = randn(p,q);

        %construct a_{t,R} based on mol
        A = template(mol,N); 


        cov = sigma^2.*eye(N);
        
        M = 2000; % number of random examples
        p_0 = 0.5; % prior prob for noise (no signal)
        [y,tl_class] = randdata(M,A,sigma,p_0); % generate y and true labels
        
        % predict labels by minimizing distance (norm)
        if (j==1)
            [pl_class,ds_cur] = detect_max(y,A,@(y,a)-d1(y,a));
        else
            [pl_class,ds_cur] = detect_max(y,A,@(y,a)d2(y,a), tao);
        end
        
        
        fp = sum(tl_class==0 & pl_class>0)/sum(tl_class==0); % false positive rate
%         fn = sum(tl_class>0 & pl_class==0)/sum(tl_class>0); % false negative rate
        
        fps(i) = fp;
%         fns(i) = fn;
        
        %             image_show = abs(sigma-sigma_error)<1e-14;
        image_show = 0;
        
        %         C = error_matrix(tl_class,pl_class,Nc,image_show); % error matrix for (t,R) pair
        C_red = error_matrix_red((tl_class),(pl_class),Nc,Nt,image_show); % reduced error matrix for (t,R) pair
        
        [g,h1,h2,h3,h4,o1,o2,o3,o4,r] = extract_C(C_red,p,Nt,image_show);
        gs(i) = g; % avg rate for true t, wrong R per t
        h1s(i) = h1; % avg fp when R=1
        h2s(i) = h2; % avg fp when R=2
        h3s(i) = h3; % avg fp when R=3
        h4s(i) = h4; % avg fp when R=4
        rs(i) = r; % avg correct (t,R)
%         o1s(j, i) = o1; % approx fn when R=1
%         o2s(j, i) = o2; % approx fn when R=2
%         o3s(j, i) = o3; % approx fn when R=3
%         o4s(j, i) = o4; % approx fn when R=4
        
        %         Ct = error_matrix(get_tr(tl_class),get_tr(pl_class),Nt,image_show); % error matrix for translation t
        Ct_red = error_matrix_red(get_tr(tl_class),get_tr(pl_class),Nt,Nt,image_show); % error matrix for t
        
        % extract error rates from Ct_red
        [a,b,c,d,e,F] = extract_Ct(Ct_red,p);
        
%         cs(j, i) = c; % tn for noise
        ds(i) = d; % tp per translation
%         as(j, i) = a; % avg fp per translation
%         bs(j, i) = b; % fn per translation
%         es(j, i) = e; % avg misclassifing rate (overlapping) per translation
%         Fs(j, i) = F; % noise rate (non-overlapping) per translation

    end
    
    dmean(j) = mean(ds);
    bmean(j) = mean(bs);
    fpmean(j) = mean(fps);
    h1mean(j) = mean(h1s);
    h2mean(j) = mean(h2s);
    h3mean(j) = mean(h3s);
    h4mean(j) = mean(h4s);
    rmean(j) = mean(rs);
    gmean(j) = mean(gs);

end



figure;
hold on;
mol_bar = zeros(2,4);
mol_bar(:,1) = rmean;
mol_bar(:,2) = dmean-rmean;
mol_bar(:,3) = 1-dmean-bmean;
mol_bar(:,4) = bmean;
bar(1:2,mol_bar,'stacked');
names = strings(2);
for i=1:2
    names(i) = sprintf('method %d', i);
end
xticks(1:2);
xticklabels(names);
title('actual mol');
legend('correct (t,R)', 'correct t, wrong R','wrong t', 'fn');

figure;
hold on;
noise_bar = zeros(2,2);
noise_bar(:,1) = 1-fpmean;
noise_bar(:,2) = fpmean;
bar(1:2,noise_bar,'stacked');
xticks(1:2);
xticklabels(names);
title('actual noise');
legend('tn','fp');

figure;
hold on;
rot_bar = zeros(2,2);
rot_bar(:,1) = h1mean+h3mean;
rot_bar(:,2) = h2mean+h4mean;
rot_bar = rot_bar./sum(rot_bar,2); % normalize over the 4 fps (sum to 1)
bar(1:2,rot_bar,'stacked');
xticks(1:2);
xticklabels(names);
legend('fat', 'tall');
title('False Positive rate per orientation')

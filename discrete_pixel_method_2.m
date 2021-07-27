clear
addpath('utils','dp');

% PS you'll find it easier to make user flags to control which figs
% produced, or, better: make separate functions to do plots of a certain
% type, reading the generated data arrays.

seeds = 0:10;
sigma = 1;

N = 24; %N grids = # pixels in 1D
ps = [4,8];%[10,8,7,4,9];
qs = [4,2];%[2,2,3,4,3];
taos = 0.2:0.01:0.8; % 0 or 1 to decide normalize the mol or not
% norms = [1];%[1,2,3]     % three diff ways to normalize: l1,l2,linfty.

k=1;
for j=1:length(ps) % iterate thru different ratios

    p = ps(j); %molecule length (p>=q)
    q = qs(j); %1 %molecule width
    Nc = (N-p+1)*4; % number of configurations (signal)
    Nt = N-p+1; % number of translations
    

    for i=1:length(taos) % iterate thru different norms (normalize or not)
        
        
        fps = zeros(1,length(seeds));
        fns = zeros(1,length(seeds));
        tps = zeros(1,length(seeds));
        ds = zeros(1,length(seeds));
        
        for l=1:length(seeds)
            
            seed = seeds(l); % rand seed generator
            rng(seed);
            
            mol = molecule(p,q,seed); %random molecule in 2D
            tao = taos(i); % threshold
            
            
            %construct a_{t,R} based on mol
            A = template(mol,N);
            
            
            cov = sigma^2.*eye(N);
            
            M = 2000; % number of random examples
            p_0 = 0.5; % prior prob for noise (no signal)
            [y,tl_class] = randdata(M,A,sigma,p_0); % generate y and true labels
            
            %             if (abs(sigma-sigma_show)<1e-14)
            %                 plot_data_sig(y,A,tl_class); % plot y with clean signals
            %             end
            
            % predict labels by maximizing <ahat, yhat> with threshold
            pl_class = detect_max(y,A,@(y,a)d2(y,a),tao);
            
            fp = sum(tl_class==0 & pl_class>0)/sum(tl_class==0); % false positive rate
            fn = sum(tl_class>0 & pl_class==0)/sum(tl_class>0); % false negative rate
            tp = sum(tl_class>0 & pl_class==tl_class)/sum(tl_class>0); % correct rate
            
            fps(l) = fp;
            fns(l) = fn;
            tps(l) = tp;
            
            %             image_show = abs(sigma-sigma_error)<1e-14;
            image_show = 0;
            
            %         C = error_matrix(tl_class,pl_class,Nc,image_show); % error matrix for (t,R) pair
            %         C_red = error_matrix_red((tl_class),(pl_class),Nc,Nt,image_show); % reduced error matrix for (t,R) pair
            
            %         [g,h1,h2,h3,h4,o1,o2,o3,o4,r] = extract_C(C_red,p,Nt,image_show);
            %         gs(k) = g; % avg rate for true t, wrong R per t
            %         h1s(k) = h1; % avg fp when R=1
            %         h2s(k) = h2; % avg fp when R=2
            %         h3s(k) = h3; % avg fp when R=3
            %         h4s(k) = h4; % avg fp when R=4
            %         o1s(k) = o1; % approx fn when R=1
            %         o2s(k) = o2; % approx fn when R=2
            %         o3s(k) = o3; % approx fn when R=3
            %         o4s(k) = o4; % approx fn when R=4
            %         rs(j,i) = r; % avg correct (t,R)
            
            %         Ct = error_matrix(get_tr(tl_class),get_tr(pl_class),Nt,image_show); % error matrix for translation t
            Ct_red = error_matrix_red(get_tr(tl_class),get_tr(pl_class),Nt,Nt,image_show); % error matrix for t
            
            % extract error rates from Ct_red
            [a,b,c,d,e,F] = extract_Ct(Ct_red,p);
            %
            %         cs(k) = c; % tn for noise
            ds(l) = d; % tp per translation
%             as(j,i) = a; % avg fp per translation
%             bs(j,i) = b; % fn per translation
            
        end
        
        dmean(j,i) = mean(ds);
        fpmean(j,i) = mean(fps);
        fnmean(j,i) = mean(fns);
        tpmean(j,i) = mean(tps);
        
        if ((tao==0.52 && j==1) || (tao==0.51 && j==2))
            fp_show(j) = fpmean(j,i);
            fn_show(j) = fnmean(j,i);
            tp_show(j) = tpmean(j,i);
            d_show(j) = dmean(j,i);
        end
        
%         es(k) = e; % avg misclassifing rate (overlapping) per translation
%         Fs(k) = F; % noise rate (non-overlapping) per translation
% 
%         a_min(k) = min(sqrt(sum(A.^2,2))); % store nearest signal to the origin
%         k = k+1;
         
    end

end

figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),1,j);
    hold on;
    plot(taos,fpmean(j,:), '.', 'Markersize', 10);
    plot(taos,fnmean(j,:), '.', 'Markersize', 10);
    xlabel('tao');
    title(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
end
legend('actual fp','actual fn');

figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),1,j);
    hold on;
    plot(fpmean(j,:), 1-fnmean(j,:),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
    xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
    title('ROC (general fp, general tp)');
    plot(fp_show(j), 1-fn_show(j),'s', 'Markersize', 15) % desired tao
end

figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),1,j);
    hold on;
    plot(fpmean(j,:), tpmean(j,:),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
    xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
    title('ROC (general fp, correct (t,R))');
    plot(fp_show(j), tp_show(j),'s', 'Markersize', 15) % desired tao
end


figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),1,j);
    hold on;
    plot(fpmean(j,:), dmean(j,:),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
    xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
    title('ROC (general fp, correct t)');
    plot(fp_show(j), d_show(j),'s', 'Markersize', 15) % desired tao
end


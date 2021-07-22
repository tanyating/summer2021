clear
addpath('utils','dp');

% PS you'll find it easier to make user flags to control which figs
% produced, or, better: make separate functions to do plots of a certain
% type, reading the generated data arrays.

seed = 0; % rand seed generator
rng(seed);

sigma = 1;
% sigmas = 0.01:0.01:2; %0.1:0.1:10; %noise levels
% sigma_show = 0.3; % noise level to plot data
% sigma_error =1.0; % noise level to show error matrix

N = 24; %N grids = # pixels in 1D
ps = [4,8];%[10,8,7,4,9];
qs = [4,2];%[2,2,3,4,3];
taos = 0.5:0.01:0.7; % 0 or 1 to decide normalize the mol or not
% norms = [1];%[1,2,3]     % three diff ways to normalize: l1,l2,linfty.

k=1;
for j=1:length(ps) % iterate thru different ratios

    p = ps(j); %molecule length (p>=q)
    q = qs(j); %1 %molecule width
    Nc = (N-p+1)*4; % number of configurations (signal)
    Nt = N-p+1; % number of translations
    

    for i=1:length(taos) % iterate thru different norms (normalize or not)
        
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
        
        fps(j,i) = fp;
        fns(j,i) = fn;
        
        %             image_show = abs(sigma-sigma_error)<1e-14;
        image_show = 0;
        
        %         C = error_matrix(tl_class,pl_class,Nc,image_show); % error matrix for (t,R) pair
        C_red = error_matrix_red((tl_class),(pl_class),Nc,Nt,image_show); % reduced error matrix for (t,R) pair
        
%         [g,h1,h2,h3,h4,o1,o2,o3,o4] = extract_C(C_red,p,Nt,image_show);
%         gs(k) = g; % avg rate for true t, wrong R per t
%         h1s(k) = h1; % avg fp when R=1
%         h2s(k) = h2; % avg fp when R=2
%         h3s(k) = h3; % avg fp when R=3
%         h4s(k) = h4; % avg fp when R=4
%         o1s(k) = o1; % approx fn when R=1
%         o2s(k) = o2; % approx fn when R=2
%         o3s(k) = o3; % approx fn when R=3
%         o4s(k) = o4; % approx fn when R=4
        
        %         Ct = error_matrix(get_tr(tl_class),get_tr(pl_class),Nt,image_show); % error matrix for translation t
        Ct_red = error_matrix_red(get_tr(tl_class),get_tr(pl_class),Nt,Nt,image_show); % error matrix for t
        
        % extract error rates from Ct_red
%         [a,b,c,d,e,F] = extract_Ct(Ct_red,p);
%         
%         cs(k) = c; % tn for noise
%         ds(k) = d; % tp per translation
%         as(k) = a; % avg fp per translation
%         bs(k) = b; % fn per translation
%         es(k) = e; % avg misclassifing rate (overlapping) per translation
%         Fs(k) = F; % noise rate (non-overlapping) per translation
% 
%         a_min(k) = min(sqrt(sum(A.^2,2))); % store nearest signal to the origin
        k = k+1;
         
    end

end

figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),1,j);
    hold on;
    plot(taos,fps(j,:), '.', 'Markersize', 10);
    plot(taos,fns(j,:), '.', 'Markersize', 10);
    xlabel('tao');
    title(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
end
legend('actual fp','actual fn');


% (sub)plots
% tt = min(sigmas):0.01:max(sigmas);
% 
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(taos)
%         ax = subplot(length(ps),length(taos),k);
%         hold on;
%         b_near = a_min(k)/2; % nearest decision boundary for x.a / ||a||
%         % plot fp
% %         plot(ax, tt, 1/2-1/2*erf(b_near./(sqrt(2).*tt)), '--', 'Linewidth', 2);hold on;
% %         plot(ax, tt, min((1/2-1/2*erf(b_near./(sqrt(2).*tt)))*Nc,1), '--', 'Linewidth', 2);
%         plot(ax, sigmas, fps(k,:), '.', 'Markersize', 10);
%         
%         % plot fn
%         plot(ax, sigmas, fns(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');
% %         legend('LB for fp','UB for fp','actual fp','actual fn');
%         vline(1,'k:','sigma=1');
%         
%         title(sprintf('p=%d q=%d tao=%.2f', ps(j),qs(j),taos(i)));
%         k=k+1;
%     end
% end
% legend('actual fp','actual fn');


% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(taos)
%         ax = subplot(length(ps),length(taos),k);
%         hold on;
%         plot(ax, sigmas, cs(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, ds(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');
%         
%         title(sprintf('p=%d q=%d tao=%.2f', ps(j),qs(j),taos(i)));
%         k=k+1;
%     end    
% end
% sgtitle('true rates');
% legend('c: tn for noise','d: tp per t for mol');
% 
% 
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(taos)
%         ax = subplot(length(ps),length(taos),k);
%         hold on;
%         plot(ax, sigmas, as(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, h1s(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, h2s(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, h3s(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, h4s(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');
%         
%         title(sprintf('p=%d q=%d tao=%.2f', ps(j),qs(j),taos(i)));
%         k=k+1;
%     end
% end
% sgtitle('false positive rates (molecule vs. noise)');
% legend('a: avg fp per t (noise)', 'h1: avg fp per t when R=1', 'h2: avg fp per t when R=2', 'h3: avg fp per t when R=3', 'h4: avg fp per t when R=4');
% 
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(taos)
%         ax = subplot(length(ps),length(taos),k);
%         hold on;
%         plot(ax, sigmas, bs(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, o1s(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, o2s(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, o3s(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, o4s(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');
%         
%         title(sprintf('p=%d q=%d tao=%.2f', ps(j),qs(j),taos(i)));
%         k=k+1;
%     end
% end
% sgtitle('false negative rates (molecule vs. noise)');
% legend('b: fn per t (noise)', 'o1: fn per t when R=1', 'o2: fn per t when R=2', 'o3: fn per t when R=3', 'o4: fn per t when R=4');
% 
% 
% 
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(taos)
%         ax = subplot(length(ps),length(taos),k);
%         hold on;
%         plot(ax, sigmas, es(k,:), '.', 'Markersize', 10);
%         hold on;
%         plot(ax, sigmas, Fs(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, gs(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');
%         
%         title(sprintf('p=%d q=%d tao=%.2f', ps(j),qs(j),taos(i)));
%         k=k+1;
%     end
% end
% sgtitle('mis-classifying rates (mol vs. mol)');
% legend('e: avg overlapping mis-class rate per t', 'F: general non-overlapping mis-class rate per t', 'g: avg true t, wrong R rate per t');

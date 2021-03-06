clear
addpath('utils','dp');

seed = 0;   % rand seed generator
rng(seed);

sigmas = 5; %0.1:0.1:10; %noise levels
sigma_show = 0.3; % noise level to plot data
sigma_error =1.0; % noise level to show error matrix

N = 96; %N grids = # pixels in 1D
ps = [16,32];%[10,8,7,4,9];
qs = [16,8];%[2,2,3,4,3];
norms = [1]; % 0 or 1 to decide normalize the mol or not
% norms = [1];%[1,2,3]     % three diff ways to normalize: l1,l2,linfty.

k=1;
for j=1:length(ps) % iterate thru different ratios

    p = ps(j); %molecule length (p>=q)
    q = qs(j); %1 %molecule width
    Nc = (N-p+1)*4; % number of configurations (signal)
    Nt = N-p+1; % number of translations
    

%     for i=1:length(norms) % iterate thru different norms (normalize or not)
        
        mol = molecule(p,q,seed); %random molecule in 2D
%         mol = randn(p,q);

        %construct a_{t,R} based on mol
        A = template(mol,N); 

        for l=1:length(sigmas)
            sigma = sigmas(l);
            cov = sigma^2.*eye(N);

            M = 5000; % number of random examples
            p_0 = 0.5; % prior prob for noise (no signal)
            [y,tl_class] = randdata(M,A,sigma,p_0,0); % generate y and true labels
            tl_pairs = inverse_map(tl_class);

%             if (abs(sigma-sigma_show)<1e-14) 
%                 plot_data_sig(y,A,tl_class); % plot y with clean signals
%             end

            % predict labels by minimizing distance (norm)
            pl_class = detect_min(y,A,@(y,a)d1(y,a));
            pl_pairs = inverse_map(pl_class);
            
            % plot true t vs. detected t (location)
            figure;
            plot(tl_pairs(tl_class~=0 & pl_class~=0 ,1), pl_pairs(tl_class~=0 & pl_class~=0,1),'r.', 'Markersize', 10);
            title(sprintf('p=%d q=%d', ps(j),qs(j)));
            xlabel('true t');
            ylabel('detected t');
            axis tight;

            % FPRs
            fps(k, l) = sum(tl_class==0 & pl_class>0)/sum(tl_class==0); % false positive rate
            fat_fps(k,l) = sum(tl_class==0 & (pl_pairs(:,2)==1 | pl_pairs(:,2)==3))/sum(tl_class==0); % fat fp rate
            tall_fps(k,l) = sum(tl_class==0 & (pl_pairs(:,2)==2 | pl_pairs(:,2)==4))/sum(tl_class==0); % tall fp rate
%             fns(k, l) = sum(tl_class>0 & pl_class==0)/sum(tl_class>0); % false negative rate

            % TPRs
            tp1s(k,l) = sum(tl_class>0 & pl_class>0)/sum(tl_class>0); % overall tp rate
            fat_tp1s(k,l) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & pl_class>0)/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat tp rate
            tall_tp1s(k,l) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & pl_class>0)/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall tp rate
            
            % Correct (t,R) rates
            tp2s(k,l) = sum(tl_class>0 & (pl_class==tl_class))/sum(tl_class>0); % correct (t,R) rate
            fat_tp2s(k,l) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & (pl_class==tl_class))/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat correct (t,R)
            tall_tp2s(k,l) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & (pl_class==tl_class))/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall correct (t,R)
            
            % Correct t rates
            tp3s(k,l) = sum(tl_class>0 & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_class>0); % correct t rate
            fat_tp3s(k,l) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat correct (t,R)
            tall_tp3s(k,l) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall correct (t,R)

            
            % extract error matrix C
%             image_show = abs(sigma-sigma_error)<1e-14;
            image_show = 0;

            %         C = error_matrix(tl_class,pl_class,Nc,image_show); % error matrix for (t,R) pair
            C_red = error_matrix_red((tl_class),(pl_class),Nc,Nt,image_show); % reduced error matrix for (t,R) pair

            [g,h1,h2,h3,h4,o1,o2,o3,o4] = extract_C(C_red,p,Nt,image_show);
            gs(k,l) = g; % avg rate for true t, wrong R per t
            h1s(k,l) = h1; % avg fp when R=1
            h2s(k,l) = h2; % avg fp when R=2
            h3s(k,l) = h3; % avg fp when R=3
            h4s(k,l) = h4; % avg fp when R=4
            o1s(k,l) = o1; % approx fn when R=1
            o2s(k,l) = o2; % approx fn when R=2
            o3s(k,l) = o3; % approx fn when R=3
            o4s(k,l) = o4; % approx fn when R=4

            
            % extract trans-wise error matrix Ct
            %         Ct = error_matrix(get_tr(tl_class),get_tr(pl_class),Nt,image_show); % error matrix for translation t
            Ct_red = error_matrix_red(get_tr(tl_class),get_tr(pl_class),Nt,Nt,image_show); % error matrix for t

            % extract error rates from Ct_red
            [a,b,c,d,e,F] = extract_Ct(Ct_red,p);

            cs(k,l) = c; % tn for noise
            ds(k,l) = d; % tp per translation
            as(k,l) = a; % avg fp per translation
            bs(k,l) = b; % fn per translation
            es(k,l) = e; % avg misclassifing rate (overlapping) per translation
            Fs(k,l) = F; % noise rate (non-overlapping) per translation

        end
        
%     a_min(k) = min(sqrt(sum(A.^2,2))); % store nearest signal to the origin   
    k = k+1;
         
%     end

end

% (sub)plots
tt = min(sigmas):0.01:max(sigmas);

%% plot FP & FN
figure;
k=1;
for j=1:length(ps)
    for i=1:length(norms)
        ax = subplot(length(ps),length(norms),k);
        hold on;
        % plot fp
%         plot(ax, sigmas, fps(k,:), '.', 'Markersize', 10);
        
        % plot fn
%         plot(ax, sigmas, 1-tp1s(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');

        bar(1,fps(k,:));
        bar(2,1-tp1s(k,:));
        xticks(1:2);
        xticklabels({'actual fp','actual fn'});

%         legend('LB for fp','UB for fp','actual fp','actual fn');
        % vline(sigma_show,'k:','sigma shown in fig.1');
        
        title(sprintf('p=%d q=%d norm=%d', ps(j),qs(j),norms(i)));
        k=k+1;
    end
end
% legend('actual fp','actual fn');


%% plot TP & TN
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(norms)
%         ax = subplot(length(ps),length(norms),k);
%         hold on;
%         plot(ax, sigmas, cs(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, ds(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');
%         
%         title(sprintf('p=%d q=%d norm=%d', ps(j),qs(j),norms(i)));
%         k=k+1;
%     end    
% end
% sgtitle('true rates');
% legend('c: tn for noise','d: tp per t for mol');

%% plot FP per orientation
figure;
k=1;
for j=1:length(ps)
    for i=1:length(norms)
        ax = subplot(length(ps),length(norms),k);
        hold on;
%         plot(ax, sigmas, fps(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, fat_fps(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, tall_fps(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');

        bar(1,fps(k,:));
        bar(2,fat_fps(k,:));
        bar(3,tall_fps(k,:));
        xticks(1:3);
        xticklabels({'overall', 'fat', 'tall'});
        
%         title(sprintf('p=%d q=%d norm=%d', ps(j),qs(j),norms(i)));
        k=k+1;
    end
end
sgtitle('False positive rates');
% legend('overall', 'fat', 'tall');

%% plot FP per t per orientation
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(norms)
%         ax = subplot(length(ps),length(norms),k);
%         hold on;
%         plot(ax, sigmas, as(k,:), '.', 'Markersize', 10);
%         plot(ax, sigmas, (h1s(k,:)+h3s(k,:)), '.', 'Markersize', 10);
%         plot(ax, sigmas, (h2s(k,:)+h4s(k,:)), '.', 'Markersize', 10);
% %         plot(ax, sigmas, h3s(k,:), '.', 'Markersize', 10);
% %         plot(ax, sigmas, h4s(k,:), '.', 'Markersize', 10);
%         xlabel('\sigma');
%         
%         title(sprintf('p=%d q=%d norm=%d', ps(j),qs(j),norms(i)));
%         k=k+1;
%     end
% end
% sgtitle('avg false positive rates per t (molecule vs. noise)');
% legend('a: overall', 'h1+h3: fat', 'h2+h4: tall');

%% plot FN per t per orientation
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(norms)
%         ax = subplot(length(ps),length(norms),k);
%         hold on;
% %         plot(ax, sigmas, bs(k,:), '.', 'Markersize', 10);
% %         plot(ax, sigmas, (o1s(k,:)+o3s(k,:))./2, '.', 'Markersize', 10);
% %         plot(ax, sigmas, (o2s(k,:)+o4s(k,:))./2, '.', 'Markersize', 10);
% 
%         bar(bs(k,:));
%         bar((o1s(k,:)+o3s(k,:))./2);
%         bar((o2s(k,:)+o4s(k,:))./2);
% %         xlabel('\sigma');
%         
%         title(sprintf('p=%d q=%d norm=%d', ps(j),qs(j),norms(i)));
%         k=k+1;
%     end
% end
% sgtitle('False negative rates (molecule vs. noise)');
% legend('overall', 'fat', 'tall');
% % legend('b: overall', 'mean(o1,o3): fat', 'mean(o2,o4): tall');


%% plot mis-classifying rates
% figure;
% k=1;
% for j=1:length(ps)
%     for i=1:length(norms)
%         ax = subplot(length(ps),length(norms),k);
%         hold on;
% %         plot(ax, sigmas, es(k,:), '.', 'Markersize', 10);
% %         plot(ax, sigmas, Fs(k,:), '.', 'Markersize', 10);
% %         plot(ax, sigmas, gs(k,:), '.', 'Markersize', 10);
% %         xlabel('\sigma');
%         
%         bar(es(k,:));
%         bar(Fs(k,:));
%         bar(gs(k,:));
% 
%         title(sprintf('p=%d q=%d norm=%d', ps(j),qs(j),norms(i)));
%         k=k+1;
%     end
% end
% sgtitle('mis-classifying rates (mol vs. mol)');
% legend('e: avg overlapping mis-class rate per t', 'F: general non-overlapping mis-class rate per t', 'g: avg true t, wrong R rate per t');


%% plot ROC curves
figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),3,k);
    hold on;
    plot(fps(j,end), tp1s(j,end),'.', 'Markersize', 10);
    plot(fat_fps(j,end), fat_tp1s(j,end),'.', 'Markersize', 10);
    plot(tall_fps(j,end), tall_tp1s(j,end),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
%     xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigmas(end)));
    xlabel('general fp');
    ylabel('general tp');
    
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(fps(j,end), tp3s(j,end),'.', 'Markersize', 10);
    plot(fat_fps(j,end), fat_tp3s(j,end),'.', 'Markersize', 10);
    plot(tall_fps(j,end), tall_tp3s(j,end),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
%     xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigmas(end)));
    xlabel('general fp');
    ylabel('correct t');
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(fps(j,end), tp2s(j,end),'.', 'Markersize', 10);
    plot(fat_fps(j,end), fat_tp2s(j,end),'.', 'Markersize', 10);
    plot(tall_fps(j,end), tall_tp2s(j,end),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
%     xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigmas(end)));
    xlabel('general fp');
    ylabel('correct (t,R)');
    k=k+1;
    
end
legend('overall','fat','tall');
sgtitle('ROC');
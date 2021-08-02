clear
addpath('utils','dp');

% PS you'll find it easier to make user flags to control which figs
% produced, or, better: make separate functions to do plots of a certain
% type, reading the generated data arrays.

seed = 0;
rng(seed);

sigma = 2;

N = 48; %N grids = # pixels in 1D
ps = [8,16];%[10,8,7,4,9];
qs = [8,4];%[2,2,3,4,3];
taos = 0.2:0.01:0.8; % 0 or 1 to decide normalize the mol or not
% norms = [1];%[1,2,3]     % three diff ways to normalize: l1,l2,linfty.

tao_picks = [0.44 0.44]; % optimal threshold for overall fp vs. fn
i_picks = (tao_picks-taos(1))./0.01 + 1;
i_picks = int16(i_picks);

fat_tao_picks = [0.4 0.4]; % optimal threshold for fat fp vs. fn
fat_i_picks = (fat_tao_picks-taos(1))./0.01 + 1;
fat_i_picks = int16(fat_i_picks);

tall_tao_picks = [0.49 0.44]; % optimal threshold for tall fp vs. fn
tall_i_picks = (tall_tao_picks-taos(1))./0.01 + 1;
tall_i_picks = int16(tall_i_picks);

k=1;
for j=1:length(ps) % iterate thru different ratios

    p = ps(j); %molecule length (p>=q)
    q = qs(j); %1 %molecule width
    Nc = (N-p+1)*4; % number of configurations (signal)
    Nt = N-p+1; % number of translations
    
    mol = molecule(p,q,seed); %random molecule in 2D
    
    %construct a_{t,R} based on mol
    A = template(mol,N);
    
    
    cov = sigma^2.*eye(N);
    
    M = 2000; % number of random examples
    p_0 = 0.5; % prior prob for noise (no signal)
    [y,tl_class] = randdata(M,A,sigma,p_0); % generate y and true labels
    tl_pairs = inverse_map(tl_class);
    

    for i=1:length(taos) % iterate thru different norms (normalize or not)
        
        tao = taos(i); % threshold
        
        % predict labels by maximizing <ahat, yhat> with threshold
        pl_class = detect_max(y,A,@(y,a)d2(y,a),tao);
        pl_pairs = inverse_map(pl_class);
        
        fps(j,i) = sum(tl_class==0 & pl_class>0)/sum(tl_class==0); % fp rate
        fat_fps(j,i) = sum(tl_class==0 & (pl_pairs(:,2)==1 | pl_pairs(:,2)==3))/sum(tl_class==0); % fat fp rate
        tall_fps(j,i) = sum(tl_class==0 & (pl_pairs(:,2)==2 | pl_pairs(:,2)==4))/sum(tl_class==0); % tall fp rate
        
        
        tp1s(j,i) = sum(tl_class>0 & pl_class>0)/sum(tl_class>0); % overall tp rate
        fat_tp1s(j,i) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & pl_class>0)/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat tp rate
        tall_tp1s(j,i) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & pl_class>0)/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall tp rate
        
        tp2s(j,i) = sum(tl_class>0 & (pl_class==tl_class))/sum(tl_class>0); % correct (t,R) rate
        fat_tp2s(j,i) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & (pl_class==tl_class))/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat correct (t,R)
        tall_tp2s(j,i) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & (pl_class==tl_class))/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall correct (t,R)
        
        tp3s(j,i) = sum(tl_class>0 & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_class>0); % correct t rate
        fat_tp3s(j,i) = sum((tl_pairs(:,2)==1 | tl_pairs(:,2)==3) & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_pairs(:,2)==1 | tl_pairs(:,2)==3); % fat correct (t,R)
        tall_tp3s(j,i) = sum((tl_pairs(:,2)==2 | tl_pairs(:,2)==4) & (pl_pairs(:,1)==tl_pairs(:,1)))/sum(tl_pairs(:,2)==2 | tl_pairs(:,2)==4); % tall correct (t,R)
        


    end

end

figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),3,k);
    hold on;
    plot(taos,fps(j,:), '.', 'Markersize', 10);
    plot(taos,fat_fps(j,:), '.', 'Markersize', 10);
    plot(taos,tall_fps(j,:), '.', 'Markersize', 10);
    plot(taos,1-tp1s(j,:), '.', 'Markersize', 10);
    xlabel(sprintf('tao (p=%d q=%d sigma=%.2f)', ps(j),qs(j),sigma));
    title('Overall mol');
    legend('overall fp','fat fp','tall fp','overall fn');
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(taos,fat_fps(j,:), '.', 'Markersize', 10);
    plot(taos,1-fat_tp1s(j,:), '.', 'Markersize', 10);
    xlabel(sprintf('tao (p=%d q=%d sigma=%.2f)', ps(j),qs(j),sigma));
    title('Fat mol');
    legend('fat fp','fat fn');
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(taos,tall_fps(j,:), '.', 'Markersize', 10);
    plot(taos,1-tall_tp1s(j,:), '.', 'Markersize', 10);
    xlabel(sprintf('tao (p=%d q=%d sigma=%.2f)', ps(j),qs(j),sigma));
    title('Tall mol');
    legend('tall fp','tall fn');
    k=k+1;
end
% legend('overall fp','fat fp','tall fp','overall fn')
% legend('overall fp','fat fp','tall fp','overall fn','fat fn','tall fn');

figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),3,k);
    hold on;
    plot(taos,tp1s(j,:), '.', 'Markersize', 10);
    plot(taos,fat_tp1s(j,:), '.', 'Markersize', 10);
    plot(taos,tall_tp1s(j,:), '.', 'Markersize', 10);
    xlabel(sprintf('tao (p=%d q=%d sigma=%.2f)', ps(j),qs(j),sigma));
    title('general TP rates');
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(taos,tp2s(j,:), '.', 'Markersize', 10);
    plot(taos,fat_tp2s(j,:), '.', 'Markersize', 10);
    plot(taos,tall_tp2s(j,:), '.', 'Markersize', 10);
    xlabel(sprintf('tao (p=%d q=%d sigma=%.2f)', ps(j),qs(j),sigma));
    title('correct (t,R) rates');
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(taos,tp3s(j,:), '.', 'Markersize', 10);
    plot(taos,fat_tp3s(j,:), '.', 'Markersize', 10);
    plot(taos,tall_tp3s(j,:), '.', 'Markersize', 10);
    xlabel(sprintf('tao (p=%d q=%d sigma=%.2f)', ps(j),qs(j),sigma));
    title('correct t rates');
    k=k+1;
end
sgtitle('TP rates');
legend('overall','fat','tall');

figure;
k=1;
for j=1:length(ps)
    ax = subplot(length(ps),3,k);
    hold on;
    plot(fps(j,:), tp1s(j,:),'.', 'Markersize', 10);
    plot(fat_fps(j,:), fat_tp1s(j,:),'.', 'Markersize', 10);
    plot(tall_fps(j,:), tall_tp1s(j,:),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
    xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
    title('ROC (general fp, general tp)');
    
    % desired tao
    plot(fps(j,i_picks(j)), tp1s(j,i_picks(j)),'s', 'Markersize', 15);
    plot(fat_fps(j,fat_i_picks(j)), fat_tp1s(j,fat_i_picks(j)),'s', 'Markersize', 15);
    plot(tall_fps(j,tall_i_picks(j)), tall_tp1s(j,tall_i_picks(j)),'s', 'Markersize', 15);
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(fps(j,:), tp2s(j,:),'.', 'Markersize', 10);
    plot(fat_fps(j,:), fat_tp2s(j,:),'.', 'Markersize', 10);
    plot(tall_fps(j,:), tall_tp2s(j,:),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
    xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
    title('ROC (general fp, correct (t,R))');
    
    % desired tao
    plot(fps(j,i_picks(j)), tp2s(j,i_picks(j)),'s', 'Markersize', 15);
    plot(fat_fps(j,fat_i_picks(j)), fat_tp2s(j,fat_i_picks(j)),'s', 'Markersize', 15);
    plot(tall_fps(j,tall_i_picks(j)), tall_tp2s(j,tall_i_picks(j)),'s', 'Markersize', 15);
    k=k+1;
    
    ax = subplot(length(ps),3,k);
    hold on;
    plot(fps(j,:), tp3s(j,:),'.', 'Markersize', 10);
    plot(fat_fps(j,:), fat_tp3s(j,:),'.', 'Markersize', 10);
    plot(tall_fps(j,:), tall_tp3s(j,:),'.', 'Markersize', 10);
    xx=0:0.01:1;
    plot(xx,xx);
    xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
    title('ROC (general fp, correct t)');
    
    % desired tao
    plot(fps(j,i_picks(j)), tp3s(j,i_picks(j)),'s', 'Markersize', 15);
    plot(fat_fps(j,fat_i_picks(j)), fat_tp3s(j,fat_i_picks(j)),'s', 'Markersize', 15);
    plot(tall_fps(j,tall_i_picks(j)), tall_tp3s(j,tall_i_picks(j)),'s', 'Markersize', 15);
    k=k+1;
end
legend('overall','fat','tall');

for j=1:length(ps)
    fprintf('When p=%d q=%d sigma=%.2f:\n', ps(j),qs(j),sigma);
    fprintf('optimal threshold for overall molecule: tao=%.2f\n', taos(i_picks(j)));
    fprintf('optimal threshold for fat molecule: tao=%.2f\n', taos(fat_i_picks(j)));
    fprintf('optimal threshold for tall molecule: tao=%.2f\n', taos(tall_i_picks(j)));
end

% figure;
% 
% for j=1:length(ps)
%     ax = subplot(length(ps),1,j);
%     hold on;
%     plot(fps(j,:), tp2s(j,:),'.', 'Markersize', 10);
%     plot(fat_fps(j,:), fat_tp2s(j,:),'.', 'Markersize', 10);
%     plot(tall_fps(j,:), tall_tp2s(j,:),'.', 'Markersize', 10);
%     xx=0:0.01:1;
%     plot(xx,xx);
%     xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
%     title('ROC (general fp, correct (t,R))');
%     
%     % desired tao
%     plot(fps(j,i_picks(j)), tp2s(j,i_picks(j)),'s', 'Markersize', 15);
%     plot(fat_fps(j,fat_i_picks(j)), fat_tp2s(j,fat_i_picks(j)),'s', 'Markersize', 15);
%     plot(tall_fps(j,tall_i_picks(j)), tall_tp2s(j,tall_i_picks(j)),'s', 'Markersize', 15);
% end
% legend('overall','fat','tall');
% 
% figure;
% 
% for j=1:length(ps)
%     ax = subplot(length(ps),1,j);
%     hold on;
%     plot(fps(j,:), tp3s(j,:),'.', 'Markersize', 10);
%     plot(fat_fps(j,:), fat_tp3s(j,:),'.', 'Markersize', 10);
%     plot(tall_fps(j,:), tall_tp3s(j,:),'.', 'Markersize', 10);
%     xx=0:0.01:1;
%     plot(xx,xx);
%     xlabel(sprintf('p=%d q=%d sigma=%.2f', ps(j),qs(j),sigma));
%     title('ROC (general fp, correct t)');
%     
%     % desired tao
%     plot(fps(j,i_picks(j)), tp3s(j,i_picks(j)),'s', 'Markersize', 15);
%     plot(fat_fps(j,fat_i_picks(j)), fat_tp3s(j,fat_i_picks(j)),'s', 'Markersize', 15);
%     plot(tall_fps(j,tall_i_picks(j)), tall_tp3s(j,tall_i_picks(j)),'s', 'Markersize', 15);
% end
% legend('overall','fat','tall');
% 

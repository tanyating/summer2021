clear
addpath('utils','dp');

N = 10; %N grids = # pixels in 1D
p = 2; %molecule length (p>=q)
q = 2; %1 %molecule width
Nc = (N-p+1)*4; % number of configurations (signal)
Nt = N-p+1; % number of translations
mol = molecule(p,q) %.*randi(10,q,p); %random molecule in 2D

%construct a_{t,R} based on mol
A = template(mol,N); 

sigmas = 0.01:0.01:2; %0.1:0.1:10; %noise
fps = zeros(size(sigmas));
fns = zeros(size(sigmas));

sigma_show = 0.3; % noise level to plot data
sigma_error =1.0; % noise level to show error matrix


for l=1:length(sigmas)
    sigma = sigmas(l);
    cov = sigma^2.*eye(N);
    
    M = 10000; % number of random examples
    p_0 = 0.5; % prior prob for noise (no signal)
    [y,tl_class] = randdata(M,A,sigma,p_0); % generate y and true labels
    % *** suggest don't use tl_pair, instead extract from tl_class when need
    if (abs(sigma-sigma_show)<1e-14) 
        plot_data_sig(y,A,tl_class); % plot y with clean signals
    end

    % predict labels by minimizing distance (norm)
    tao = 0.7;
    pl_class = detect_max(y,A,@(y,a)d2(y,a),tao);

    C = error_matrix(tl_class,pl_class,Nc,0);
    fp = sum(C(1,2:end)); % false positive rate
    fn = sum(tl_class>0 & pl_class==0)/sum(tl_class>0); % false negative rate
    % fn = sum(C(2:end,1))/sum(C(2:end,:),'all')
    
    fps(l) = fp;
    fns(l) = fn;
    if (abs(sigma-sigma_error)<1e-14) % visualize error matrix
        C = error_matrix(tl_class,pl_class,Nc); % error matrix for (t,R) pair
        C_red = error_matrix_red((tl_class),(pl_class),Nc,Nt)
        
        Ct = error_matrix(get_tr(tl_class),get_tr(pl_class),Nt) % error matrix for translation t
        Ct_red = error_matrix_red(get_tr(tl_class),get_tr(pl_class),Nt,Nt)
    end
    
end

tt = min(sigmas):0.01:max(sigmas);

%plot fp
figure;
plot(sigmas, fps, '.', 'Markersize', 10);

%plot fn
hold on;
plot(sigmas, fns, '.', 'Markersize', 10);
xlabel('\sigma');
legend('actual fp','actual fn');
% vline(sigma_show,'k:','sigma shown in fig.1');
hold off;

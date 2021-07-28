clear;
addpath('utils','dp');

raw = 0;   % AHB added: 0: use best dist as in meth2,   1: use all dists (Manas) giving tail as in Rickgauer Fig 1e.
sigma = 1; %level of noise
seeds = 0:10;

N = 24; %N grids = # pixels in 1D
p = 8; %molecule length (p>=q)
q = 2; %1 %molecule width
Nc = (N-p+1)*4; % number of configurations (signal)
Nt = N-p+1; % number of translations

for i=1:length(seeds)

    seed =  seeds(i); rng(seed);
    
    mol = molecule(p,q,seed); %random molecule in 2D
    
    %construct a_{t,R} based on mol
    A = template(mol,N);
    
    M = 10000; %number of noise vectors
    cov = sigma^2.*eye(N);
    p_0 = 0.5; % prior prob for noise (no signal)
    [x,tl_class] = randdata(M,A,sigma,1); % generate noise
    [y,tl_class] = randdata(M,A,sigma,0); % generate signals
    
    
    [x_pl,x_d2] =  detect_max(x,A,@(y,a)d2(y,a));
    [y_pl,y_d2] =  detect_max(y,A,@(y,a)d2(y,a));
    
    x_d2s(i,:) = transpose(x_d2(:));
    y_d2s(i,:) = transpose(y_d2(:));
    
    max_x_d2s(i,:) = max(x_d2,[],2);
    max_y_d2s(i,:) = max(y_d2,[],2);

end

nbins=100;



figure;
if raw
  histogram(x_d2s(:),nbins); hold on;
  histogram(y_d2s(:),nbins);
  title('all inner products');
else
  histogram(mean(max_x_d2s,1),nbins); hold on;
  histogram(mean(max_y_d2s,1),nbins);
  title('max inner product');
end
legend('no signal', 'signal');
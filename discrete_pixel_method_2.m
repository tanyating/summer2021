addpath('utils');
clf;

Rs = [0, pi/2, pi, pi*3/2]; %4 rotations
N = 10; %N grids
p = 2; %molecule length (p>q)
q = 1; %molecule width
mol = rand(q,p).*randi(10,q,p); %random molecules in 2D
ts = 0:N-p; %N-p+1 translations
A = zeros(N-p+1,4,N); %store all a_{t,R} (index:translation, rotation, vector)

%construct a_{t,R} based on mol
for i=0:N-p
    for j=1:4
        tmp=sum(rot90(mol,j-1),1);
        A(1+i,j,1+i:i+length(tmp)) = tmp;
    end
end

sigmas = 0.1:0.01:1; %noise
fps = zeros(size(sigmas));
fns = zeros(size(sigmas));

sigma_show = 0.1;

for l=1:length(sigmas)
    sigma = sigmas(l);
    cov = sigma^2.*eye(N);

    M = 10000; %number of random examples
    y = mvnrnd(zeros(N,1), cov, M); %add noise to each example
    truelabel = zeros(M,1);
    filter = rand(M,1); %uniform distribution 
    p_0 = 0.5; %prior for no signal
    tmp = (1-p_0)/(4*(N-p+1)); %equally likely for all classes with signal
%     tmp = 1/(4*(N-p+1)+1);
    truelabel(filter>p_0) = 1;
%     truelabel(filter>tmp) = 1;

     % plot no signal class (guassian blobs) in 2D
    if (N==2 && sigma==sigma_show)
        plot(y((truelabel==0),1),y((truelabel==0),2),'.','Markersize',10);
        hold on;
    end

    %assign each random vector with signal based on true label
    a = zeros(1,N);
    k=0;
    for i=0:N-p
        for j=1:4
            a(1,:) = A(i+1,j,:); %a_{t,R}
            y(filter>(p_0+k*tmp) & filter<=(p_0+(k+1)*tmp), :) = y(filter>(p_0+k*tmp) & filter<=(p_0+(k+1)*tmp), :) + a;
            
            % plot signal classes (guassian blobs) in 2D
            if (N==2 && sigma==sigma_show)
                plot(y(filter>(p_0+k*tmp) & filter<=(p_0+(k+1)*tmp),1),y(filter>(p_0+k*tmp) & filter<=(p_0+(k+1)*tmp),2),'.','Markersize',10);
            end
            
            k = k+1;
        end
    end


    % predict labels by minimizing distance (norm)
    predlabel = zeros(M,1);

    % set up a threshold based on prior (bias)
    tao = 0.8;
    yhat = y-mean(y,2);
    yhat = yhat./sqrt(sum(yhat.^2,2));
    
    a = zeros(1,N);
    ahat = zeros(1,N);
    cur_max = yhat*transpose(ahat);
    for i=0:N-p
       cur_prod = zeros(M,4);
       for j=1:4
           a(1,:) = A(i+1,j,:);
           ahat = (a-mean(a))/norm(a-mean(a));
           cur_prod(:,j) = yhat*transpose(ahat); %inner product <xhat,a_{t,R}hat>
       end
       % threshold predictor
       predlabel(max(cur_prod,[],2)>tao) = 1; %f(t)>tao, then candidate
       
       % maximize <xhat,a_{t,R}hat>
%        maxprod_vec = max(cur_prod,[],2);
%        predlabel(maxprod_vec>cur_max) = 1;
%        cur_max(maxprod_vec>cur_max) = maxprod_vec(maxprod_vec>cur_max);
    end
    
    
    fp = sum((truelabel-predlabel)==-1)/sum(truelabel==0);
    fn = sum((truelabel-predlabel)==1)/sum(truelabel==1);
    
    fps(l) = fp;
    fns(l) = fn;


end

% fprintf('Max dot product results in translation %d and rotation %.4f\n', max_i, Rs(max_j));


tt = min(sigmas):0.01:max(sigmas);
% b = norm(a)/2+tt.^2.*log(p_0/p_a)./norm(a); % decision boundary for x.a / ||a||

%plot fp
figure(2);
% plot(tt, 1/2-1/2*erf(b./(sqrt(2).*tt)), '-', 'Linewidth', 2);
hold on;
plot(sigmas, fps, '.', 'Markersize', 10);

%plot fn
% figure(3);
hold on;
% plot(tt, 1/2-1/2*erf(-(b-norm(a))./(sqrt(2).*tt)), '-', 'Linewidth', 2);
plot(sigmas, fns, '.', 'Markersize', 10);
xlabel('\sigma');
% legend('expexted fp', 'actual fp', 'expexted fn', 'actual fn');
legend('actual fp','actual fn');
% vline(sigma_show,'k:','sigma shown in fig.1');
hold off;

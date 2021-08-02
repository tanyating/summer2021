% script to test linear discriminant analysis (LDA) for presence/absence of a
% signal vector in R^d, given additive iid Gaussian noise, and compare failure
% rates to analytic formulae, over various noise levels. Maximum a-posteriori
% (MAP) selection is used over the two labels (absent=0, present=1),
% where posterior \propto likelihood * prior. The decision hyperplane is shown
% for one noise level.

% Notes: As opposed to ML (max likelihood) which has decision plane as the
% perp bisector of 0 and a (the signal vector), MAP can at high noise levels
% give a boundary that is way outside the region between 0 and a. This is as
% expected since it is telling you that the likelihood is uninformative, and
% the prior dominates the MAP choice. Basically if pi_1 > 1/2, for noise->infty
% LDA always predicts no signal, so FPR->0 and FNR->1.

% to fix: whether class labels called 0,1 or 1,2 throughout.

% Tanya Wang. June 2021. Input from Alex Barnett and Manas Rachh.

addpath('utils')   % include local stuff incl stats/ML toolbox copies

a = [5,3]; %signal   (note it's at a generic angle, tests more than (5,5) does)
sigmas = 1:0.5:20; %noise
sigma_show = 1;    % must be exactly in the above list
d = length(a); %dimension
p_0 = 0.5; %prior belief in no signal = its incoming rate = pi_1 in notes
p_a = 1-p_0; %prior belief in signal = also its incoming rate = pi_2 in notes

fps = zeros(size(sigmas));
fns = zeros(size(sigmas));

for i=1:length(sigmas) %noise
    
    sigma = sigmas(i);
    cov = sigma^2*eye(d);

    % N random examples with prior
    N = 10000;
    truelabel = rand(N,1)>p_0; % label 0 or 1
    x = mvnrnd(zeros(d,1), cov, N);
    x(truelabel,:) = x(truelabel,:) + a;
    
    % plot gaussian blobs in 2-D
    if (d==2 && sigma==sigma_show)
        figure(1);clf
        plot(x(truelabel==0,1),x(truelabel==0,2), 'r.', 'Markersize', 5);
        hold on;
        plot(x(truelabel==1,1),x(truelabel==1,2), 'b.', 'Markersize', 5);
        xx = -1:0.1:1;
        hold on;
        plot(a(1)/2+sigma^2*log(p_0/p_a)*cos(atan(a(2)/a(1)))-a(2).*xx, a(2)/2+sigma^2*log(p_0/p_a)*sin(atan(a(2)/a(1)))+a(1).*xx, 'k--', 'Linewidth', 2, 'Markersize', 10);
        plot(0,0,'r*'); plot(a(1),a(2),'b*');   % add the means (true signals)
        legend('class 1', 'class 2', 'decision boundary');
        axis equal         % to get aspect ratio right
    end
    
    predlabel = zeros(N,1);
    
    % (1) max likelihood predictor
    % predlabel(mvnpdf(x, zeros(1,d), cov)*p_0 < mvnpdf(x, a, cov)*p_a) = 1;

    % (2) decision boundary predictor
    predlabel(x*transpose(a) > norm(a)^2/2 + sigma^2*log(p_0/p_a)) = 1;

    fp = sum((truelabel-predlabel)==-1)/sum(truelabel==0); %false positive rate
    fn = sum((truelabel-predlabel)==1)/sum(truelabel==1); %false negative rate

    %fprintf('false positive rate is %.6f\n', fp);
    %fprintf('false negative rate is %.6f\n', fn);
    
    
    fps(i) = fp;
    fns(i) = fn;
  
end



tt = min(sigmas):0.01:max(sigmas);
b = norm(a)/2+tt.^2.*log(p_0/p_a)./norm(a); % decision boundary for x.a / ||a||

%plot fp
figure(2);
plot(tt, 1/2-1/2*erf(b./(sqrt(2).*tt)), '-', 'Linewidth', 2);
hold on;
plot(sigmas, fps, '.', 'Markersize', 10);
% xlabel('\sigma');
% legend('expexted fp', 'actual fp');

%plot fn
% figure(3);
hold on;
plot(tt, 1/2-1/2*erf(-(b-norm(a))./(sqrt(2).*tt)), '-', 'Linewidth', 2);
plot(sigmas, fns, '.', 'Markersize', 10);
xlabel('\sigma');
legend('expexted fp', 'actual fp', 'expexted fn', 'actual fn');
vline(sigma_show,'k:','sigma shown in fig.1');
hold off;



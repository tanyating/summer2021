a = [5,5]; %signal
sigmas = 1:0.5:50; %noise
d = length(a); %dimension
p_0 = 0.5; %prior belief in no signal
p_a = 1-p_0; %prior belief in signal

fps = zeros(size(sigmas));
fns = zeros(size(sigmas));

for i=1:length(sigmas) %noise
    
    sigma = sigmas(i);
    cov = sigma^2*eye(d);

    % N random examples with prior
    N = 10000;
    truelabel = rand(N,1)>p_0; % label 0 or 1
    x = mvnrnd(zeros(d,1), cov, N);
    x(truelabel==1,:) = x(truelabel==1,:) + a;
    
    % plot gaussian blobs in 2-D
    if (d==2 && sigma==1)
        plot(x(truelabel==0,1),x(truelabel==0,2), 'r.', 'Markersize', 5);
        hold on;
        plot(x(truelabel==1,1),x(truelabel==1,2), 'b.', 'Markersize', 5);
        xx = -1:0.1:1;
        hold on;
        plot(a(1)/2+sigma^2*log(p_0/p_a)*cos(atan(a(2)/a(1)))-a(2).*xx, a(2)/2+sigma^2*log(p_0/p_a)*sin(atan(a(2)/a(1)))+a(1).*xx, 'k--', 'Linewidth', 2, 'Markersize', 10);
        legend('class 1', 'class 2', 'decision boundary');
    end
    
    predlabel = zeros(N,1);
    
    % (1) max likelihood predictor
    % predlabel(mvnpdf(x, zeros(1,d), cov)*p_0 < mvnpdf(x, a, cov)*p_a) = 1;

    % (2) decision boundary predictor
    predlabel(x*transpose(a) > norm(a)^2/2 + sigma^2*log(p_0/p_a)) = 1;

    fp = sum(truelabel-predlabel==-1)/sum(truelabel==0); %false positive rate
    fn = sum(truelabel-predlabel==1)/sum(truelabel==1); %false negative rate

    %fprintf('false positive rate is %.6f\n', fp);
    %fprintf('false negative rate is %.6f\n', fn);
    
    
    fps(i) = fp;
    fns(i) = fn;
  
end



tt = min(sigmas):0.01:max(sigmas);
b = norm(a)/2+tt.^2.*log(p_0/p_a)./norm(a); %decision boundary

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
hold off;



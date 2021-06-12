a = [5,5]; %signal
sigmas = 1:0.5:50;
%sigmas = [0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]; %noise
d = length(a); %dimension
p_0 = 0.8; %prior belief in no signal
p_a = 1-p_0; %prior belief in signal

fps = zeros(size(sigmas));
fns = zeros(size(sigmas));

for i=1:length(sigmas) %noise
    
    sigma = sigmas(i);
    cov = sigma^2*eye(d);

    % random examples with no signal
    N = 10000;
    N1 = int16(N*p_0);
    N2 = N-N1;
    x1 = zeros(1,d)+mvnrnd(zeros(d,1), cov, N1);
    labels1_1 = zeros(N1,1); %label 0 indicates no signal, label 1 indicates signal

    % random examples with no signal
    x2 = a+mvnrnd(zeros(d,1), cov, N2);
    labels2_1 = zeros(N2,1);
    
    % plot gaussian blobs in 2-D
    if (d==2 && sigma==1)
        plot(x1(:,1),x1(:,2), 'r.', 'Markersize', 5);
        hold on;
        plot(x2(:,1),x2(:,2), 'b.', 'Markersize', 5);
        xx = -1:0.1:1;
        hold on;
        plot(a(1)/2+sigma^2*log(p_0/p_a)*cos(atan(a(2)/a(1)))-a(2).*xx, a(2)/2+sigma^2*log(p_0/p_a)*sin(atan(a(2)/a(1)))+a(1).*xx, 'k--', 'Linewidth', 2, 'Markersize', 10);
        legend('class 1', 'class 2', 'decision boundary');
    end

    % max likelihood predictor
    labels1_1(mvnpdf(x1, zeros(1,d), cov)*p_0 < mvnpdf(x1, a, cov)*p_a) = 1;
    labels2_1(mvnpdf(x2, zeros(1,d), cov)*p_0 < mvnpdf(x2, a, cov)*p_a) = 1;

    % for i=1:N
    %     if (mvnpdf(x1, zeros(d,1), cov)*p_0 >= mvnpdf(x1, a, cov)*p_a)
    %          labels(i)=1;
    %         fprintf("no signal\n");
    %     else
    %          labels(i)=2;
    %         fprintf("signal\n");
    %     end
    % end

    labels1_2 = zeros(N1,1);
    labels2_2 = zeros(N2,1);

    % decision boundary predictor
    labels1_2(x1*transpose(a) > norm(a)^2/2 + sigma^2*log(p_0/p_a)) = 1;
    labels2_2(x2*transpose(a) > norm(a)^2/2 + sigma^2*log(p_0/p_a)) = 1;

    % for i=1:N
    %     if (transpose(x1(i))*a <= norm(a)^2/2 + sigma^2*log(p_0/p_a))
    %         labels2(i)=1;
    %         fprintf("no signal\n");
    %     else
    %         labels2(i)=2;
    %         fprintf("signal\n");
    %     end   
    % end

    fp = sum(labels1_1)/length(labels1_1); %false positive rate
    fn = 1-sum(labels2_1)/length(labels2_1); %false negative rate

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



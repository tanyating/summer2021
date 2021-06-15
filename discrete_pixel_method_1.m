Rs = [0, pi/2, pi, pi*3/2]; %4 rotations
N = 10;
m = [1 2 0; 0 10 0];
[q,p] = size(m); %molecule length and width
ts = 0:N-p; %N-p+1 translations
A = zeros(N-p+1,4,N); %store all a_{t,R} (index:translation, rotation, vector)

%construct a_{t,R} based on m
for i=0:N-p
    for j=1:4
        tmp=sum(rot90(m,j-1),1);
        A(1+i,j,1+i:i+length(tmp)) = tmp;
    end
end

sigma = 0.1;

a = zeros(1,N);
a(1,:) = A(5+1,1,:);
y = mvnrnd(a,sigma^2*eye(N),1); %given vector with random noise sigma

% (1) minimize |y-a_{t,R}|
min_i = 0;
min_j = 1;
a = zeros(1,length(A(1,1,:)));
a(1,:) = A(1,1,:);
cur_min= norm(y-a);
for i=0:N-p
    for j=1:4
        a = zeros(1,length(A(i+1,j,:)));
        a(1,:) = A(i+1,j,:);
        if (norm(y-a) < cur_min)
            min_i = i;
            min_j = j;
            cur_min= norm(y-a);
        end
    end
end

fprintf('Min norm search results in translation %d and rotation %.4f\n', min_i, Rs(min_j));

% (2) maximize likelihood (multivariate Gaussian)
max_i = 0;
max_j = 1;
a = zeros(1,length(A(1,1,:)));
a(1,:) = A(i,1,:);
cur_max = mvnpdf(y,a,sigma^2*eye(N));
for i=0:N-p
    for j=1:4
        a = zeros(1,length(A(i+1,j,:)));
        a(1,:) = A(i+1,j,:);
        if (mvnpdf(y,a,sigma^2*eye(N)) > cur_max)
            max_i = i;
            max_j = j;
            cur_max = mvnpdf(y,a,sigma^2*eye(N));
        end
    end
end

fprintf('Max likelihood search results in translation %d and rotation %.4f\n', max_i, Rs(max_j));

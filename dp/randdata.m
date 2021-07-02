function [y,tl_class] = randdata(M,A,sigma,p0)
% RANDDATA  Generate random data matrix of size M*N (M instances of dim N)
% given template A of size (N-p+1)*4*N, prior p0, and noise level sigma.
%
% [y,tl_class] = randdata(M,A,sigma,p0) returns data y of size M*N, 
% with approximately p0*M instances of noise drawn from Guassion N(0,sigma^2 I), 
% (1-p0)/Nc * M instances drawn from each configuration N(a_{t,R}, sigma^2 I), 
%
% and the tl_class returns M*1 array indicating the true class 
% (unique combination of translation i and rotation j indexed by
% k=i*4+j in {0,..,Nc}) for each instance.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_randdata; return; end
if nargin<4, p0=0.5; end

Nc = size(A,1); % number of configurations
N = size(A,2); % N grids
cov = sigma^2.*eye(N);

y = mvnrnd(zeros(N,1), cov, M); % add noise to each example
tl_class = zeros(M,1); % true class

filter = rand(M,1); %uniform distribution
tmp = (1-p0)/Nc; %equally likely for each configuration

%assign each random vector with signal based on true label
a = zeros(1,N);
for k=1:Nc
        a(1,:) = A(k,:); % a_{t,R}
        
        lb = filter>(p0+(k-1)*tmp); % probability lower bound
        ub = filter<=(p0+(k)*tmp); % probability upper bound
        y(lb & ub, :) = y(lb & ub, :) + a; % assign current signal
        tl_class(lb & ub) = k; % assign unique class
            
end



%%%%%%%%
function test_randdata 

N = 3;
p = 2;
Nc = (N-p+1)*4;
A = rand(Nc,N); % random template
p0 = 0.5;
sigma = 0.1;
M = 100; % generates 100 random samples
[y,tl_class] = randdata(M,A,sigma,p0);
% check if the pair labels and unique labels match
tl_pair = inverse_map(tl_class);
if (max(tl_pair(:,1))>N-p)
    error('incorrect label');
elseif (abs(sum(tl_class==0)-50)>10)
    error('incorrect instances of noise');
else
    'ok'
end



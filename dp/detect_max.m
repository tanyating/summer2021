function [pl_class,pl_pair] = detect_max(y,A,d)
% DETECT_MAX  Detect/predict class for each row of data y (M*N) based on
% template A ((N-p+1)*4*N) by maximizing "distance" (function) d.
%
% pl_class =  detect_max(y,A,d) returns the predicted labels pl_class of M*1 array indicating
% detected class (unique combination of translation i and rotation j) for each
% instance.
%
% [pl_class,pl_pair] = detect_max(y,A,d) ... also returns M*2 array pl_pair
% of true pairs of labels, each row [i,j] indicating translation i and jth
% rotation detected for each instance.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_detect_max; return; end
if (size(A,3)~=size(y,2)) error('N dimension must be the same for A and y!'); end

M = size(y,1); % number of configurations
N = size(A,3); % N grids
p = N+1-size(A,1); % size(A,1) = N-p+1

pl_class = zeros(M,1); % initialize predicted class
pl_pair = zeros(M,2); % initialize predicted pairs

a = zeros(1,N); % start from "distance" with the origin (no signal)
cur_max = d(y,a);
k=1;
for i=0:N-p
    for j=1:4
        a(1,:) = A(i+1,j,:); % "distance" to current a_{t,R}
        tmp = d(y,a);
        pl_class(tmp>cur_max) = k; % predict with max "distance" label
        pl_pair(tmp>cur_max, 1) = i;
        pl_pair(tmp>cur_max, 2) = j;
        cur_max(tmp>cur_max) = tmp(tmp>cur_max); % swap current max "distance"
        k=k+1;
    end
end


%%%%%%%%
function test_detect_max 
A = zeros(1,4,3); % simple template with 4 configurations
A(1,1,:) = [1 1 1];
A(1,2,:) = [2 2 2];
A(1,3,:) = [3 3 3];
A(1,4,:) = [4 4 4];
y = [0 0 0;1 1 1]; % signals without noise
pl_class =  detect_max(y,A,@(y,a)(-d1(y,a))); % min norm
if (pl_class(1) ~= 0 | pl_class(2) ~= 1)
    error('failed');
end

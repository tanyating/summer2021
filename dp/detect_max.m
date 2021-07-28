function [pl_class,ds] = detect_max(y,A,d,tao)
% DETECT_MAX  Detect/predict class for each row of data y (M*N) based on
% template A ((N-p+1)*4*N) by maximizing "distance" (function) d. If there
% is a threshold tao, the predicted label is positive if the max-distance >
% tao and returns the max-distance set of translation and rotation, else
% returns no signal (noise).
%
% pl_class =  detect_max(y,A,d) returns the predicted labels pl_class of M*1 array indicating
% detected class (unique combination of translation i and rotation j) for each
% instance.
%
% [pl_class,ds] = detect_max(y,A,d) ... also returns M*Nc array of
% "distances", each entry of ith row and kth column indicating
% d(y_i,a_k), which is the "distance" between the ith instance and kth
% template.
%
% [pl_class,ds] = detect_max(y,A,d,tao) ... also checks whether the
% max-distance > tao to determine positive (signal) or not.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_detect_max; return; end
if (size(A,2)~=size(y,2)) error('N dimension must be the same for A and y!'); end

M = size(y,1); % number of instances
N = size(A,2); % N grids
Nc = size(A,1); % number of configurations

pl_class = zeros(M,1); % initialize predicted class
ds = zeros(M,Nc); % initialize "distances"

a = zeros(1,N); % start from "distance" with the origin (no signal)
cur_max = d(y,a);

for k=1:Nc
    
        a(1,:) = A(k,:); % current a_{t,R}
        tmp = d(y,a); % "distance" to current a_{t,R}
        ds(:,k) = tmp;
        pl_class(tmp>cur_max) = k; % predict with max "distance" label
        cur_max(tmp>cur_max) = tmp(tmp>cur_max); % swap current max "distance"
    
end

if (nargin>3) % checks with threshold
        pl_class(cur_max<tao) = 0;
end


%%%%%%%%
function test_detect_max 
A = [1 1 1; 2 2 2; 3 3 3; 4 4 4]; % simple template with 4 configurations
y = [0 0 0;1 1 1]; % signals without noise
pl_class =  detect_max(y,A,@(y,a)(-d1(y,a))); % min norm
if (pl_class(1) ~= 0 | pl_class(2) ~= 1)
    error('failed');
else
    'ok'
end

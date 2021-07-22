function [g,h1,h2,h3,h4,o1,o2,o3,o4] = extract_C(C_red,p,Nt,center_show)
% EXTRACT-CT  Extracts different error rates per translation-rotation pair
% from the reduced confusion matrix C given the maximum length of the
% molecule p, and the total number of translations Nt=N-p+1.
%
% [g,h1,h2,h3,h4] = extract_C(C_red,p,Nt) returns g as the average rate for true t wrong R per translation,
% h1 as average false positive rate when predicted R=1,
% h2 as average false positive rate when predicted R=2,
% h3 as average false positive rate when predicted R=3,
% h4 as average false positive rate when predicted R=4,
% o1 as approximate false negative rate when predicted R=1,
% o2 as approximate false negative rate when predicted R=2,
% o3 as approximate false negative rate when predicted R=3,
% o4 as approximate false negative rate when predicted R=4.
%
% If the flag center_show is set to 1 (default to be 0), the function also plots the error
% matrix for the the central translation with 4 rotations (overlapped
% center).
%
% Without arguments, a self-test is done.
%

% Tanya 7/12/21.

if nargin==0, test_extract_C; return; end
if nargin<4, center_show=0; end

if (size(C_red,1)~=5) error('Incorrect size for reduced error matrix'); end

% get the overlapped region
overlap_width = p-1;
overlap_center = map_class([ceil(Nt/2)-1 1])+1; % overlapped central t w/ 1st R
left = overlap_center-overlap_width*4;
right = overlap_center+3+overlap_width*4;
overlap = C_red(2:end,left:right);

% visualize overlapped central t w/ 4 R's (4-by-4 matrix)
if (center_show)
    figure;imagesc(overlap(:,overlap_width*4+1:overlap_width*4+4));title('center of overlap');colorbar; colormap(jet(256));
    xlabel('pred label'); ylabel('true label');axis equal;caxis([0 1]);
end

ovtruetrans = overlap(:,overlap_width*4+1:overlap_width*4+4);   % avoids 'all'
g = (sum(ovtruetrans(:))-sum(diag(ovtruetrans)))/(16-4); % avg rate for true t, wrong R per t
h = mean(reshape(C_red(1,2:end),4,[]),2); % vector of avg fp per R (length 4 array)
h1 = h(1); % avg fp when R=1
h2 = h(2); % avg fp when R=2
h3 = h(3); % avg fp when R=3
h4 = h(4); % avg fp when R=4
o1 = C_red(2,1); % approx fn when R=1
o2 = C_red(3,1); % approx fn when R=1
o3 = C_red(4,1); % approx fn when R=1
o4 = C_red(5,1); % approx fn when R=1


%%%%%%%%
function test_extract_C
Nt = 5;
C_red = zeros(5, Nt*4+1);
p = 2;
C_red(:,1) = 1;
[g,h1,h2,h3,h4] = extract_C(C_red,p,Nt);
if (g~=0 | h1~=0 | h2~=0 | h3~=0 | h4~=0)
    error('failed');
else
    'ok'
end

Nt=5;
C_red = zeros(5, Nt*4+1);
p = 2;
C_red(:,1:5) = eye(5);
[g,h1,h2,h3,h4] = extract_C(C_red,p,Nt);
if (g~=0 | h1~=0 | h2~=0 | h3~=0 | h4~=0)
    error('failed');
else
    'ok'
end
% 
% Ct_red = [1 0 0; 1 0 0];
% p = 1;
% [a,b,c,d,e,F] = extract_Ct(Ct_red,p)
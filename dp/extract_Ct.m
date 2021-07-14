function [a,b,c,d,e,F] = extract_Ct(Ct_red,p)
% EXTRACT-CT  Extracts different error rates per translation from the
% reduced translation-wise confusion matrix  Ct given the maximum length of the
% molecule p.
%
% [a,b,c,d,e,F] = extract_Ct(Ct_red,p) returns a as the average false positive rate per translation,
% b as the rough false negative rate per translation,
% c as the true negative rate for noise,
% d as the rough true positive rate per translation,
% e as the average mis-classifying rate per overlapping/touching translation,
% f as the average noise (mis-classifying) rate per non-overlapping translation.
%
% Without arguments, a self-test is done.
%

% Tanya 7/12/21.

if nargin==0, test_extract_Ct; return; end

Nt = size(Ct_red,2)-1;

if (size(Ct_red,1)~=2) error('Incorrect size for reduced translation error matrix'); end

% overlapped region t-wise
overlap_width = p-1;
overlap_center = ceil(Nt/2)+1; % central t
left = overlap_center-overlap_width;
right = overlap_center+overlap_width;
overlap = Ct_red(2,left:right);

c = Ct_red(1,1); % tn for noise
d = overlap(overlap_width+1); % tp per translation
a = mean(Ct_red(1,2:end)); % avg fp per translation
b = Ct_red(2,1); % fn per translation
e = mean(cat(2,overlap(1:overlap_width),overlap(overlap_width+2:end))); % avg misclassifing rate (overlapping) per translation
F = mean(cat(2,Ct_red(2,2:left-1),Ct_red(2,right+1:end))); % noise rate (non-overlapping) per translation


%%%%%%%%
function test_extract_Ct
Ct_red = [1 0 0 0 0 0; 0 0 0 1 0 0];
p = 2;
[a,b,c,d,e,F] = extract_Ct(Ct_red,p);
if (a~=0 | b~=0 | c~=1 | d~=1 | e~=0 | F~=0)
    error('failed');
else
    'ok'
end

Ct_red = [0 0.1 0.1 0.2 0.3 0.3; 0.1 0.1 0.2 0.3 0.2 0.1];
p = 2;
[a,b,c,d,e,F] = extract_Ct(Ct_red,p);
if (a~=0.2 | b~=0.1 | c~=0 | d~=0.3 | e~=0.2 | F~=0.1)
    error('failed');
else
    'ok'
end
% 
% Ct_red = [1 0 0; 1 0 0];
% p = 1;
% [a,b,c,d,e,F] = extract_Ct(Ct_red,p)
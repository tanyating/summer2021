function [A,AA] = template(mol,N)
% TEMPLATE  Generate (N-p+1)*4 number of 1D templates (R^N)
% based on a molecule mol (p*q) in 2D.
%
% A = template(mol,N) generates a [(N-p+1),4,N] table A
% such that A(i,j,:) corrresponds to a pure signal with translation i and rotation j.
%
% [A,AA] = template(mol,N) ... also returns a [(N-p+1)*4,N] matrix AA
% such that each row of AA corrresponds to a pure signal.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_template; return; end

p = max(size(mol)); % length of input molecule
if (N<p) error('N must be greater or equal to either side of molecule!'); end

Nc = (N-p+1)*4; % total number of configurations
A = zeros(N-p+1,4,N); %store all a_{t,R} (index:translation, rotation, vector)

for i=0:N-p % iterate over each translation
    for j=1:4 % iterate over each rotaion
        tmp=sum(rot90(mol,j-1),1);
        A(1+i,j,1+i:i+length(tmp)) = tmp;
    end
end

AA = reshape(A, [Nc N]); % stack configuration vectors as rows


%%%%%%%%
function test_template 
[A,AA] = template(1,1); % simple template with N=1, and Nc=1*4=4
if (any(A~=[1 1 1 1]) | any(AA~=[1 1 1 1]))
    error('failed');
end
try, A = template(rand(1,4),3); catch me, ['ok: ',me.message], end
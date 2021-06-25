function d = d1(y,a)
% D1  Compute "distance 1" (method 1: norm) between data y (M*N) and signal a (1*N). 
%
% d1(y,a) calculates |y_i-a| for each ith row of y,
% and returns vector d of size M*1.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_d1; return; end

if (size(a,1)~=1) error('a must be of size 1*N!'); end
if (size(y,2)~=size(a,2)) error('a and y must have same 2nd dimension'); end

d = sqrt(sum((y-a).^2,2));

%%%%%%%%
function test_d1    % throw it a simple pair (p,q) with/out seed
try, d = d1([1,2],[0,0]); catch me, ['ok: ',me.message], end
try, d = d1([1,2],[0,0,1]); catch me, ['ok: ',me.message], end % error
try, d = d1([1,2;3,4],[0,0]); catch me, ['ok: ',me.message], end
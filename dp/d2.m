function d = d2(y,a)
% D2  Compute "distance 2" (method 2: <ahat,yhat>) between data y (M*N) and signal a (1*N). 
%
% d2(y,a) calculates dot product of yhat and ahat for each ith row of y,
% and returns vector d of size M*1.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_d2; return; end

if (size(a,1)~=1) error('a must be of size 1*N!'); end
if (size(y,2)~=size(a,2)) error('a and y must have same 2nd dimension'); end

yhat = y-mean(y,2);
ynorm = sqrt(sum(yhat.^2,2));
yhat = yhat(ynorm~=0,:)./ynorm(ynorm~=0,:); % normalize each row of y
ahat = a-mean(a); 
anorm = norm(ahat);
if (anorm~=0) 
    ahat = ahat/anorm; % normalize a
end

d = yhat*transpose(ahat);

%%%%%%%%
function test_d2    
d = d2([1,-1;-1,1],[1,-1]);
if (abs(d(1)-1)>1e-6 || abs(d(2)+1)>1e-6)
    error('failed');
else
    'ok'
end
try, d = d2([1,2],[0,0,1]); catch me, ['ok: ',me.message], end
try, d = d2([0,0],[1,2;3,4]); catch me, ['ok: ',me.message], end

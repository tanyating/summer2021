function mol = molecule(p,q,seed)
% MOLECULE  Generate random molecule matrix of size q*p (assume p>q, p as length, q as width) 
% with certain random generator seed.
%
% mol = molecule(p,q) generates random matrix of size p*q 
% by fixing default seed=0.
%
% mol = molecule(p,q,seed) resets the random generator seed
% and also returns random matrix of size p*q.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_molecule; return; end
if nargin<3, seed=0; end

if (p<q) error('p must be greater or equal to q!'); end

rng(seed); % fix seed

mol = rand(q,p);


%%%%%%%%
function test_molecule    % throw it a simple pair (p,q) with/out seed
try, mol = molecule(2,2); catch me, ['ok: ',me.message], end
try, mol = molecule(2,2,1); catch me, ['ok: ',me.message], end
try, mol = molecule(3,2); catch me, ['ok: ',me.message], end
try, mol = molecule(3,2,10); catch me, ['ok: ',me.message], end
try, mol = molecule(2,3); catch me, ['ok: ',me.message], end
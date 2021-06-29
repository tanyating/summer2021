function tl_pairs = inverse_map(tl_class)
% INVERSE_MAP  Maps a set of single indices of a label to combinations (rows) of translation i and rotation j.
%
% tl_pairs = inverse_map(tl_class) maps each row k to i=k/4, j=k%4
% if k=0, maps to no signal [0,0].
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_inverse_map; return; end

if (any(tl_class<0)) error('any label k must be greater or equal to 0!'); end
if (size(tl_class,2)~=1) error('any label k must be 1 dimension!'); end

tl_pairs = zeros(length(tl_class),2);
tl_pairs(:,2) = mod(tl_class,4);
tl_pairs(tl_pairs(:,2)==0,2) = 4;
tl_pairs(:,1) = (tl_class-tl_pairs(:,2))./4;

tl_pairs(tl_class==0,1) = 0;
tl_pairs(tl_class==0,2) = 0;


%%%%%%%%
function test_inverse_map    % throw it a simple pair (i,j)
tl_pairs = inverse_map([5;11]);
if (any(tl_pairs(1,:)~=[1 1]) | any(tl_pairs(2,:)~=[2 3]))
    error('failed');
end
try, tl_pairs = inverse_map([5,11]); catch me, ['ok: ',me.message], end
try, tl_pairs = inverse_map([-1;1]); catch me, ['ok: ',me.message], end 

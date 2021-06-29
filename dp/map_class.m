function tl_class = map_class(tl_pairs)
% MAP_CLASS  Maps a matrix with 2 columns of translations and rotations to single indices.
%
% k = map_class(tl_pairs) maps each row [i,j] of tl_pairs to k=i*4+j,
% if i=0 and j=0, maps to no signal k=0.
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_map; return; end

if (size(tl_pairs,2)~=2) error('there should be 2 columns for translation and rotation'); end
if (any(tl_pairs(:,2)>4)) error('rotation must be smaller or equal to 4!'); end

tl_class = tl_pairs(:,1).*4 + tl_pairs(:,2);
tl_class(tl_pairs(:,1)==0 & tl_pairs(:,2)==0) = 0;


%%%%%%%%
function test_map    % throw it a simple pair (i,j)
tl_class = map_class([1,1;2,3]);
if (tl_class(1)~=5 | tl_class(2)~=11)
    error('failed');
end
try, tl_class = map_class([1,5]); catch me, ['ok: ',me.message], end
try, tl_class = map_class([2]); catch me, ['ok: ',me.message], end
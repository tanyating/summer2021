function tr_class = get_tr(input_class)
% GET_TR  transforms an M-by-1 array of labels (unique for each pair of
% translation and rotation) into an an M-by-1 array of labels unique for
% translation only.
%
% tr_class = get_tr(input_class) first maps each label (row) k of input_class into
% corresponding pair of translation and rotation [i,j], then returns i+1 as
% a unique index for the instance's translation (returns 0 if the row is [0
% 0] as no signal).
%
% Without arguments, a self-test is done.
%

% Tanya 6/30/21.

if nargin==0, test_get_tr; return; end

input_pair = inverse_map(input_class); % get translation and rotations pair [i,j] for each label
tr_class = input_pair(:,1)+1; % get translation index [i+1]
tr_class(input_class==0) = 0; % assign 0 (no signal) class


%%%%%%%%
function test_get_tr
input_pair = [0 0; 0 1; 1 1; 2 2; 3 3];
input_class = map_class(input_pair);
tr_class = get_tr(input_class);
if (any(tr_class ~= [0;1;2;3;4]))
    error('failed');
else
    'ok'
end
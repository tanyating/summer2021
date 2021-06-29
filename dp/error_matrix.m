function C = error_matrix(tl_class,pl_class,Nc,image_show)
% ERROR_MAX  Compute confusion/error matrix based on true labels of size M*1
% and predicted labels of size M*1.
%
% C = error_matrix(tl_class,pl_class) returns the error matrix of size
% (Nc+1)*(Nc+1), each entry Cij indicating the relative occurence
% (frequency) of an instance with true label i, but detected with label j.
%
% C = error_matrix(tl_class,pl_class,Nc,image_show)... also plots the error
% matrix if image_show is true (default to be 1)
%
% Without arguments, a self-test is done.
%

% Tanya 6/29/21.

if nargin==0, test_error_matrix; return; end
if nargin<4, image_show=1; end

if (size(tl_class,1)~=size(pl_class,1)) error('Size of true labels and predicted labels should equal!'); end

M = size(tl_class,1); % number of instances (labels)

C = zeros(Nc+1,Nc+1); % Nc configurations of signal and 1 no-signal

for i=1:M % iterate over each instance
    tl_cur = tl_class(i); % current true label
    pl_cur = pl_class(i); % current predicted label
    C(tl_cur+1,pl_cur+1) = C(tl_cur+1,pl_cur+1)+1; % add to error matrix
end

C = C./sum(C,2); % normalize C along each row

if (image_show)
    figure;imagesc(C);title('Confusion matrix');colorbar; colormap(gray(256)); % visualize error matrix
end




%%%%%%%%
function test_error_matrix
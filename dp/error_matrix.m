function C = error_matrix(tl_class,pl_class,Nc,image_show)
% ERROR_MATRIX  Compute confusion/error matrix based on true labels of size M*1
% and predicted labels of size M*1.
%
% C = error_matrix(tl_class,pl_class,Nc) returns the error matrix of size
% (Nc+1)*(Nc+1), each entry Cij indicating the relative occurence
% (frequency) of an instance with true label i, but detected with label j.
% This means each row gives the fractions of classifications given to the
% true label i.
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

% C = C./M; % normalize by total number of instances M
C = C./sum(C,2); % normalize C along each row

if (image_show)
    figure;imagesc(C);title('Confusion matrix (row-normalized)');colorbar; colormap(jet(256)); % visualize error matrix
    xlabel('pred label'); ylabel('true label');axis equal;caxis([0 1]);
end


%%%%%%%%
function test_error_matrix
Nc = 5;
X = randi(Nc+1,1e3,1)-1;
C = error_matrix(X,X,Nc);
if norm(C-eye(Nc+1))==0, 'ok', else, error('failed'); end

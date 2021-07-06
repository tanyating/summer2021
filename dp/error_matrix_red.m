function C_red = error_matrix_red(tl_class,pl_class,Nc,Nt)
% ERROR_MATRIX_RED  Compute the "reduced" version of confusion/error matrix
% (only the different entries in the Toplitz matrix) based on true labels
% of size M*1 and predicted labels of size M*1.
%
% C_red = error_matrix_red(tl_class,pl_class,Nc,tr) returns the error matrix of size
% 2-by-(Nc+1) if Nc==Nt (translation-only), each entry Cij indicating the relative occurence
% (frequency) of an instance with true label i (either no signal or the
% "central" translation), and detected with label j.
%
% If Nc>Nt (translation-rotation), returns the error matrix of size
% 5-by-(Nc+1), for the rows of true label are the no signal and the
% "central" translation with 4 rotations.
%
% This means each row gives the fractions of classifications given to the
% true label i.
%
% The function also plots the "reduced" error matrix.
%
% Without arguments, a self-test is done.
%

% Tanya 6/30/21.

if nargin==0, test_error_matrix_red; return; end

if (size(tl_class,1)~=size(pl_class,1)) error('Size of true labels and predicted labels should equal!'); end

if (Nc==Nt) % translation-wise
    C_red = zeros(2,Nt+1);
    k=1; % reduced error matrix row number
    for i=[0,ceil(Nt/2)] % iterate for true label = 0 and Nt/2
        for j=0:Nt % iterate through each predicted label
            C_red(k,j+1) = sum(tl_class==i & pl_class==j); % add to error matrix
        end
        k=k+1;
    end
    
else % translation-rotation
    C_red = zeros(5,Nc+1);
    for i=0:4
        if (i==0)
            k=0; % true label = 0
        else
            k=map_class([ceil(Nt/2)-1 i]); % true label of "central" translation and rotation i
        end
        for j=0:Nc
            C_red(i+1,j+1) = sum(tl_class==k & pl_class==j); % add to error matrix
        end
    end
end


C_red = C_red./sum(C_red,2); % normalize C along each row


figure;imagesc(C_red);title('Reduced Confusion matrix (row-normalized)');colorbar; colormap(jet(256)); % visualize error matrix
xlabel('pred label'); ylabel('true label');axis equal;caxis([0 1]);



%%%%%%%%
function test_error_matrix_red
Nt = 5;
X = randi(Nt+1,1e3,1)-1;
C_red = error_matrix_red(X,X,Nt,Nt);
I = eye(Nt+1);
if norm(C_red-I([1 ceil(Nt/2)+1],:))==0, 'ok', else, error('failed'); end
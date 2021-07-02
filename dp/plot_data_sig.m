function plot_data_sig(y,A,tl_class)
% PLOT_DATA_SIG  Plot random data matrix y (M*N) and signal templates A ((N-p+)*4*N) 
% to visualize and differentiate different clusters (classes). 
%
% If N=2, plot each cluster with noiseless signal in 2D
% If N=3, plot each cluster with noiseless signal in 3D
%
% If N>3, visualize template matrix AA and data matrix y (class-based ordered/grouped).
%
% Without arguments, a self-test is done.
%

% Tanya 6/24/21.

if nargin==0, test_plot_data_sig; return; end

N = size(y,2); % dimension N
Nc = size(A,1); % number of configurations

if (N==2) % 2d plot
    figure; clf; plot(A(:,1),A(:,2),'+', 'Markersize',20); axis equal; hold on; title('clean signals'); drawnow; 
    for i=0:Nc % iterate over each tl_class
        plot(y(tl_class==i,1),y(tl_class==i,2),'.','Markersize',10);
        hold on;
    end

elseif (N==3) % 3d plot
    figure; clf; plot3(A(:,1),A(:,2),A(:,3),'+', 'Markersize',20,'DisplayName','signals'); axis equal; hold on; title('clean signals'); drawnow; axis vis3d
    for i=0:Nc % iterate over each tl_class
        plot3(y(tl_class==i,1),y(tl_class==i,2),y(tl_class==i,3),'.','Markersize',10,'DisplayName',sprintf('%d', i));
        hold on;
%         text(AA(:,1),AA(:,2),AA(:,3),sprintf('%d',i),'Fontsize',30);
    end
    legend;
    
else % higher dimension: visualize AA, y
    
    figure;imagesc(A);title('template A (pure signals)');colorbar; colormap(gray(256)); % visualize matrix AA
    
    [tl_sorted, tl_order] = sort(tl_class); % sort data y by true class (ascending order: 0,1...,Nc+1)
    y_sorted = y(tl_order,:);
    figure;imagesc(y_sorted);title('sorted random data y');colorbar; colormap(gray(256)); % visualize data clusters

end



%%%%%%%%
function test_plot_data_sig 

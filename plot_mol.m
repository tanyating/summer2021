clear;
addpath('utils','dp');

N = 96; %14; %N grids = # pixels in 1D
ps = [16,32]; %2; %molecule length (p>=q)
qs = [16,8]; %2; %1 %molecule width

% seed =0;   %  matters!

for j=1:length(ps) % different ratios
    p = ps(j);
    q = qs(j);
    for seed=[0] % different seeds
      mol = molecule(p,q,seed); %random molecule in 2D
%         figure;imagesc(mol);title('molecule');colorbar; colormap(jet(256)); % visualize molecule
        
        %construct a_{t,R} based on mol
        A = template(mol,N);
        xx = 1:N;
        
        figure;
        for z=[1 2] 
            k = map_class([(z-1)*(p+30), z]);
%             figure;imagesc(rot90(mol,z-1));title('molecule');colorbar; colormap(jet(256));caxis([-2 2]); % visualize molecule
%             stairs(xx,A(k,:),'Linewidth',2);
            bar(A(k,:));
            fprintf('When p=%d q=%d R=%d: l1 norm = %f\n', ps(j),qs(j),z, norm(A(k,:),1));
            fprintf('When p=%d q=%d R=%d: l2 norm = %f\n', ps(j),qs(j),z, norm(A(k,:),2));
            hold on;
        end
        ylim([-10 10]);
%         title(sprintf('p=%d q=%d seed=%d', ps(j),qs(j),seed));

        grid on;
%         axis equal;

    end
end
        
clear;
addpath('utils','dp');

N = 24; %14; %N grids = # pixels in 1D
ps = [4]; %2; %molecule length (p>=q)
qs = [4]; %2; %1 %molecule width

% seed =0;   %  matters!

for j=1:length(ps) % different ratios
    p = ps(j);
    q = qs(j);
    for seed=[0] % different seeds
      mol = molecule(p,q,seed); %random molecule in 2D
        figure;imagesc(mol);title('molecule');colorbar; colormap(jet(256)); % visualize molecule
        
        %construct a_{t,R} based on mol
        A = template(mol,N);
        xx = 1:N;
        
        figure;
        for z=1:4 
            k = map_class([(z-1)*(p+1), z]);
            stairs(xx,A(k,:),'Linewidth',2);
            hold on;
        end
        title(sprintf('p=%d q=%d seed=%d', ps(j),qs(j),seed));

        grid on;
        axis equal;

    end
end
        
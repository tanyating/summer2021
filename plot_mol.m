clear;
addpath('utils','dp');

N = 35; %14; %N grids = # pixels in 1D
ps = [4,8]; %2; %molecule length (p>=q)
qs = [4,2]; %2; %1 %molecule width

for j=1:length(ps) % different ratios
    p = ps(j);
    q = qs(j);
    for i=[1,0] % normalize mol or not

        mol = molecule(p,q,0,i); %random molecule in 2D
        % figure;imagesc(mol);title('molecule');colorbar; colormap(jet(256)); % visualize molecule
        
        %construct a_{t,R} based on mol
        A = template(mol,N);
        xx = 1:N;
        
        figure;
        for z=1:4 
            k = map_class([(z-1)*(p+1), z]);
            stairs(xx,A(k,:),'Linewidth',2);
            hold on;
        end
        title(sprintf('p=%d q=%d norm=%d', ps(j),qs(j),i));

        grid on;
        axis equal;

    end
end
        
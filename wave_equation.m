%clear all;
%close all;

implementations_number = 1;
Nx=500; % grid points
Nz=500;
dx = 4;
dz = 4;
L = Nx * dx;
H = Nz * dz;
src_nx = floor(Nx/2) + 1;
src_nz = floor(Nz/2) + 1;

x = linspace(0, L, Nx);
z = linspace(0, H, Nz);

% Matematicheskoe ozhidanie
mf=3000;  % middle velocity
% Standartnoe otklonenie
sf=90;
% Dlina korrelyacii
If=25;

% plot
f = figure(2);
tl = tiledlayout(f,2,2,'TileSpacing','compact');
tl.Padding = 'normal';

for iter = 1:implementations_number
    
    filename=[num2str(sf), '_', num2str(If), '_lastwaveTxx_', num2str(iter)];
    filenameExtention = [filename, '.mod'];
    fid1 = fopen(filenameExtention, 'r');
    Vp1 = fread(fid1, [Nx, Nz], 'double');
    fclose(fid1);
    
    filename=[num2str(sf), '_', num2str(If), '_lastwaveTzz_', num2str(iter)];
    filenameExtention = [filename, '.mod'];
    fid2 = fopen(filenameExtention, 'r');
    Vp2 = fread(fid2, [Nx, Nz], 'double');
    fclose(fid2);
    

    ax1 = nexttile(tl);
    imagesc(x, z, sqrt(Vp1.^2+Vp2.^2));
    colorbar
    %colormap gray
    title(['\fontsize{14}Stresses for media ', num2str(iter)])
    axis equal, axis tight
hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
    
end


xlabel(tl, 'x (m)')
ylabel(tl,'z (m)')


figure();
filename=['300_3_lastwaveTxx_main_1'];
filenameExtention = [filename, '.mod'];
fid1 = fopen(filenameExtention, 'r');
Vp1 = fread(fid1, [Nx, Nz], 'double');
fclose(fid1);
filename=['300_3_lastwaveTzz_main_1'];
filenameExtention = [filename, '.mod'];
fid1 = fopen(filenameExtention, 'r');
Vp2 = fread(fid1, [Nx, Nz], 'double');
fclose(fid1);
imagesc(x, z, sqrt(Vp1.^2+Vp2.^2));
colorbar
xlabel('x (m)')
ylabel('z (m)')
title('\fontsize{14}Stresses for initial clutter')

axis equal, axis tight

hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
lgd = legend( {'true source'});


f=figure('Color', 'black');
tl = tiledlayout(f,1,2,'TileSpacing','compact');
tl.Padding = 'normal';
filename=['clatter_main_150_2'];
filenameExtention = [filename, '.mod'];
fid1 = fopen(filenameExtention, 'r');
Vp11 = fread(fid1, [Nx, Nz], 'double');
fclose(fid1);
filename=['clatter_2_300_3'];
filenameExtention = [filename, '.mod'];
fid1 = fopen(filenameExtention, 'r');
Vp12 = fread(fid1, [Nx, Nz], 'double');
fclose(fid1);
ax1 = nexttile(tl);
summ = Vp11+Vp12;
semblance_ = zeros(Nx, Nz);
sum_max = sum(summ, 'all');
sum_max_max = max(sum(summ, 3), [], 'all');

for ii = 1:Nx
    for jj = 1:Nz
        sum_a = sum(summ(ii, jj, :));
        if (sum_a < 75 * sum_max_max / 100)
            semblance_(ii, jj) = 0;
        else
            num = sum(summ(ii, jj, :)).^2;
            den = 1 * sum(summ(ii, jj, :).^2); % ??
            semblance_(ii, jj) = num / den;
        end
    end
end
s=surf(x, z, semblance_);
colormap('cool');
camlight; shading interp;
set(s, 'facelighting', 'phong', 'facealpha', 0.7);
colorbar; 
colorbar
axis equal, axis tight
ax1 = nexttile(tl);
    filename=[num2str(sf), '_', num2str(If), '_lastwaveTxx_', num2str(iter)];
    filenameExtention = [filename, '.mod'];
    fid1 = fopen(filenameExtention, 'r');
    Vp2 = fread(fid1, [Nx, Nz], 'double');
    fclose(fid1);
    s=surf(x, z, sqrt(Vp11.^2+Vp2.^2));
    set(gca,'Color','black')
colormap('cool');
camlight; shading interp;
set(s, 'facelighting', 'phong', 'facealpha', 0.7);
colorbar; 


xlabel(tl,'x (m)')
ylabel(tl,'z (m)')
axis equal, axis tight

setColor(text, 'warn')

 

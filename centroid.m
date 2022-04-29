% close all;
% clear all;

nx = 500;
nz = 500;
dx = 4;
dz = 4;
L = nx * dx;
H = nz * dz;

nt_iter_stop = 3000; % Number of timesteps to compute
time_global = 1.5;

% Soure location
sf = 150;
If = 2;
src_nx = floor(nx/2) + 1;
src_nz = floor(nz/2) + 1;

% create grid
x = linspace(0, L, nx);
z = linspace(0, H, nz);

trms_number = 10;
trms = zeros(nx, nz, trms_number);
MAPV = zeros(nx, nz, trms_number);
PAPR = zeros(nx, nz, trms_number);
Txx = zeros(nx, nz, trms_number);
Tzz = zeros(nx, nz, trms_number);
energy = zeros(nx, nz, trms_number);
sumPAPR = zeros(nx, nz);
sumenergy = zeros(nx, nz);
sumTzz = zeros(nx, nz);
sumTxx = zeros(nx, nz);
sumtrms = zeros(nx, nz);

f = figure();
f.Position = [200 200 800 220];
tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact');
tl.Padding = 'normal';
ax1 = nexttile(tl);
for i = 1:1
    filename = [num2str(sf), '_', num2str(If), '_lastwaveTzz_het_', num2str(i)];
    filenameExtention = [filename, '.mod'];
    fidTzz = fopen(filenameExtention, 'r');
    Tzz(:, :, i) = fread(fidTzz, [nx, nz], 'double');
%     figure(i);
%     imagesc(x, z, trms(:, :, i));
    sumTzz = sumTzz + Tzz(:, :, i);
    fclose(fidTzz);
    
    filename = [num2str(sf), '_', num2str(If), '_lastwaveTxx_het_', num2str(i)];
    filenameExtention = [filename, '.mod'];
    fidTxx = fopen(filenameExtention, 'r');
    Txx(:, :, i) = fread(fidTxx, [nx, nz], 'double');
    sumTxx = sumTxx + Txx(:, :, i);
%     
    filename = [num2str(sf), '_', num2str(If), '_energy_het_', num2str(i)];
    %    filename = [num2str(sf), '_', num2str(If), '_energy_generated', num2str(i), '_20'];

    filenameExtention = [filename, '.mod'];
    fidenergy = fopen(filenameExtention, 'r');
    energy(:, :, i) = fread(fidenergy, [nx, nz], 'double');
    sumenergy = sumenergy + energy(:, :, i);
    
        filename = [num2str(sf), '_', num2str(If), '_trm_het_', num2str(i)];
    filenameExtention = [filename, '.mod'];
    fidtrms = fopen(filenameExtention, 'r');
    trms(:, :, i) = fread(fidtrms, [nx, nz], 'double');
    sumtrms = sumtrms + trms(:, :, i);
    
    
    filenamePAPR = [num2str(sf), '_', num2str(If), '_PAPR_het_', num2str(i)];
     %   filenamePAPR = [num2str(sf), '_', num2str(If), '_PAPR_generated', num2str(i), '_20'];

    filenameExtentionPAPR = [filenamePAPR, '.mod'];
    fidPAPR = fopen(filenameExtentionPAPR, 'r');
    PAPR(:, :, i) = fread(fidPAPR, [nx, nz], 'double');
    fclose(fidPAPR);
    sumPAPR = sumPAPR + PAPR(:, :, i);
end
imax =0;
jmax=0;


summ=sqrt(sumTxx.^2 + sumTzz.^2);
max_value = summ(1, 1);
for i = 1:nx
    for j = 1:nz
        if (max_value < summ(i, j))
            max_value = summ(i, j);
            imax = i;
            jmax = j;
        end
    end
end
disp(imax);
disp(jmax);
T=summ;
imagesc(x, z, T);

title('\fontsize{14}Stress components')
axis equal, axis tight

ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor='white';

hold on
plot(src_nz*dz, src_nx*dx, 'r.', ...
                    'MarkerSize', 10);
max_value = T(1, 1);
for i = 1:nx
    for j = 1:nz
        if (max_value < T(i, j))
            max_value = T(i, j);
            imax = i;
            jmax = j;
        end
    end
end
hold on
plot(jmax*dz, imax*dx,'ys', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'y');
hold on
props = regionprops(true(size(T)), T, 'WeightedCentroid');
plot(props.WeightedCentroid(1)*dz, props.WeightedCentroid(2)*dx,'gs', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'g');
ax1 = nexttile(tl);
imagesc(x, z, sumenergy);

title('\fontsize{14}Elastic energy')
axis equal, axis tight

ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor='white';
hold on
plot(src_nz*dz, src_nx*dx, 'r.', ...
                    'MarkerSize', 10);
props = regionprops(true(size(sumenergy)), sumenergy, 'WeightedCentroid');
max_value = sumenergy(1, 1);
for i = 20:nx
    for j = 1:nz
        if (max_value < sumenergy(i, j))
            max_value = sumenergy(i, j);
            imax = i;
            jmax = j;
        end
    end
end
hold on
plot(jmax*dz, imax*dx,'ys', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'y');
                hold on

plot(props.WeightedCentroid(1)*dz, props.WeightedCentroid(2)*dx,'gs', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'g');

%     
% 
% 
% 
ax1 = nexttile(tl);
imagesc(x, z, sumtrms);


title({'\fontsize{14}Time summed','elastic energy'})
axis equal, axis tight

ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor='white';
    hold on
plot(src_nz*dz, src_nx*dx, 'r.', ...
                    'MarkerSize', 10);
props = regionprops(true(size(sumtrms(30:nx,:))), sumtrms(30:nx,:), 'WeightedCentroid');                
%props = regionprops(true(size(trms)), trms, 'WeightedCentroid');
    max_value = sumtrms(1, 1);
    for i = 30:nx
        for j = 1:nz
            if (max_value < sumtrms(i, j))
                max_value = sumtrms(i, j);
                imax = i;
                jmax = j;
            end
        end
    end

hold on
plot(jmax*dz, imax*dx,'ys', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'y');

hold on
plot(props.WeightedCentroid(1)*dz, props.WeightedCentroid(2)*dx,'gs', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'g');

ax1 = nexttile(tl);
imagesc(x, z, sumPAPR);

xlabel(tl, 'x (m)')
ylabel(tl ,'z (m)')
title('\fontsize{14}PAPR')
axis equal, axis tight

ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor='white';
%sumPAPR(1:30,:,:)=0;
props = regionprops(true(size(sumPAPR(30:nx,:))), sumPAPR(30:nx,:), 'WeightedCentroid');


hold on
plot(src_nz*dz, src_nx*dx, 'r.', ...
                    'MarkerSize', 10);
max_value = sumPAPR(1, 1);
for i = 30:nx
    for j = 1:nz
        if (max_value < sumPAPR(i, j))
            max_value = sumPAPR(i, j);
            imax = i;
            jmax = j;
        end
    end
end
    hold on
plot(jmax*dz, imax*dx,'ys', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'y');

hold on
props = regionprops(true(size(sumPAPR(30:nx,:))), sumPAPR(30:nx,:), 'WeightedCentroid');
plot(props.WeightedCentroid(1)*dz, props.WeightedCentroid(2)*dx,'gs', ...
    'LineWidth', 1, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'g');
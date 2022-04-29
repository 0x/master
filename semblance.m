close all;
clear all;

nx = 500;
nz = 500;
dx = 4;
dz = 4;
L = nx * dx;
H = nz * dz;

nt_iter_stop = 3000; % Number of timesteps to compute
time_global = 1.3;

% Soure location
sf = 300;
If = 3;
src_nx = floor(nx/2) + 1;
src_nz = floor(nz/2) + 1;

% create grid
x = linspace(0, L, nx);
z = linspace(0, H, nz);

trms_number = 50;
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

f = figure(2);
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
tl.Padding = 'normal';
for i = 1:trms_number
    %     filename = [num2str(sf), '_', num2str(If), '_lastwaveTzz_', num2str(i)];
    %     filenameExtention = [filename, '.mod'];
    %     fidTzz = fopen(filenameExtention, 'r');
    %     Tzz(:, :, i) = fread(fidTzz, [nx, nz], 'double');
    %     %     figure(i);
    %     %     imagesc(x, z, trms(:, :, i));
    %     sumTzz = sumTzz + Tzz(:, :, i);
    %     fclose(fidTzz);
    %
    %     filename = [num2str(sf), '_', num2str(If), '_lastwaveTxx_', num2str(i)];
    %     filenameExtention = [filename, '.mod'];
    %     fidTxx = fopen(filenameExtention, 'r');
    %     Txx(:, :, i) = fread(fidTxx, [nx, nz], 'double');
    %     sumTxx = sumTxx + Txx(:, :, i);

    filename = [num2str(sf), '_', num2str(If), '_energy_cr_', num2str(i)];
    filenameExtention = [filename, '.mod'];
    fidenergy = fopen(filenameExtention, 'r');
    energy(:, :, i) = fread(fidenergy, [nx, nz], 'double');
    sumenergy = sumenergy + energy(:, :, i);

    %     filename = [num2str(sf), '_', num2str(If), '_trmnс_', num2str(i)];
    %     filenameExtention = [filename, '.mod'];
    %     fidtrms = fopen(filenameExtention, 'r');
    %     trms(:, :, i) = fread(fidtrms, [nx, nz], 'double');
    %     sumtrms = sumtrms + trms(:, :, i);

    filenamePAPR = [num2str(sf), '_', num2str(If), '_PAPR_cr_', num2str(i)];
    filenameExtentionPAPR = [filenamePAPR, '.mod'];
    fidPAPR = fopen(filenameExtentionPAPR, 'r');
    PAPR(:, :, i) = fread(fidPAPR, [nx, nz], 'double');
    fclose(fidPAPR);
    sumPAPR = sumPAPR + PAPR(:, :, i);
end
imax = 0;
jmax = 0;
perc = 30;


ax1 = nexttile(tl);
summ = sumTxx + sumTzz;
semblance_ = zeros(nx, nz);
sum_max = sum(summ, 'all');
sum_max_max = max(sum(summ, 3), [], 'all');
for ii = 1:nx
    for jj = 1:nz
        sum_a = sum(summ(ii, jj, :));
        if (sum_a < perc * sum_max_max / 100)
            semblance_(ii, jj) = 0;
        else
            num = sum(summ(ii, jj, :)).^2;
            den = trms_number * sum(summ(ii, jj, :).^2); % ??
            semblance_(ii, jj) = num / den;
        end
    end
end

imagesc(x, z, semblance_);
colorbar
title('\fontsize{14}Stress components')
axis equal, axis tight
grid
ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor = 'white';

hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
hold on
max_value = semblance_(1, 1);
for i = 30:nx
    for j = 1:nz
        if (max_value < semblance_(i, j))
            max_value = semblance_(i, j);
            imax = i;
            jmax = j;
        end
    end
end
plot(jmax*dz, imax*dx, 'gs', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'g', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);


ax1 = nexttile(tl);
summ = energy;
semblance_ = zeros(nx, nz);
sum_max = sum(summ, 'all');
sum_max_max = max(sum(summ, 3), [], 'all');

for ii = 1:nx
    for jj = 1:nz
        sum_a = sum(summ(ii, jj, :));
        if (sum_a < perc * sum_max_max / 100)
            semblance_(ii, jj) = 0;
        else
            num = sum(summ(ii, jj, :)).^2;
            den = trms_number * sum(summ(ii, jj, :).^2); % ??
            semblance_(ii, jj) = num / den;
        end
    end
end

imagesc(x, z, semblance_);
colorbar
title('\fontsize{14}Energy')
axis equal, axis tight
grid
ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor = 'white';
hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
hold on
max_value = semblance_(1, 1);
for i = 30:nx
    for j = 1:nz
        if (max_value < semblance_(i, j))
            max_value = semblance_(i, j);
            imax = i;
            jmax = j;
        end
    end
end
props = regionprops(true(size(semblance_)), semblance_, 'WeightedCentroid');
plot(props.WeightedCentroid(1)*dz, props.WeightedCentroid(2)*dx, 'gs', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'g', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);

ax1 = nexttile(tl);
summ = trms;
summ(1:30, :, :) = 0;
semblance_ = zeros(nx, nz);
sum_max = sum(summ, 'all');
sum_max_max = max(sum(summ, 3), [], 'all');

for ii = 1:nx
    for jj = 1:nz
        sum_a = sum(summ(ii, jj, :));
        if (sum_a < perc * sum_max_max / 100)
            semblance_(ii, jj) = 0;
        else
            num = sum(summ(ii, jj, :)).^2;
            den = trms_number * sum(summ(ii, jj, :).^2); % ??
            semblance_(ii, jj) = num / den;
        end
    end
end

imagesc(x, z, semblance_);
colorbar
title('\fontsize{14}Summed energy')
axis equal, axis tight
grid
ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor = 'white';
props = regionprops(true(size(semblance_)), semblance_, 'WeightedCentroid');
hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);

hold on
max_value = semblance_(1, 1);
for i = 30:nx
    for j = 1:nz
        if (max_value < semblance_(i, j))
            max_value = semblance_(i, j);
            imax = i;
            jmax = j;
        end
    end
end
plot(jmax*dz, imax*dx, 'gs', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'g', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);


ax1 = nexttile(tl);
summ = PAPR;
summ(1:20, :, :) = 0;
semblance_ = zeros(nx, nz);
sum_max = sum(summ, 'all');
sum_max_max = max(sum(summ, 3), [], 'all');

for ii = 1:nx
    for jj = 1:nz
        sum_a = sum(summ(ii, jj, :));
        if (sum_a < perc * sum_max_max / 100)
            semblance_(ii, jj) = 0;
        else
            num = sum(summ(ii, jj, :)).^2;
            den = trms_number * sum(summ(ii, jj, :).^2); % ??
            semblance_(ii, jj) = num / den;
        end
    end
end

imagesc(x, z, semblance_);
colorbar
xlabel(tl, 'x (m)')
ylabel(tl, 'z (m)')
title('\fontsize{14}PAPR')
axis equal, axis tight
grid
ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor = 'white';

hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);

hold on
max_value = semblance_(1, 1);
for i = 20:nx
    for j = 1:nz
        if (max_value < semblance_(i, j))
            max_value = semblance_(i, j);
            imax = i;
            jmax = j;
        end
    end
end
plot(jmax*dz, imax*dx, 'gs', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'g', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);






figure();
summ = energy;
semblance_ = zeros(nx, nz);
sum_max = sum(summ, 'all');
sum_max_max = max(sum(summ, 3), [], 'all');

for ii = 1:nx
    for jj = 1:nz
        sum_a = sum(summ(ii, jj, :));
        if (sum_a < perc * sum_max_max / 100)
            semblance_(ii, jj) = 0;
        else
            num = sum(summ(ii, jj, :)).^2;
            den = trms_number * sum(summ(ii, jj, :).^2); % ??
            semblance_(ii, jj) = num / den;
        end
    end
end
imagesc(x, z, semblance_);
colorbar
title('\fontsize{14}Energy, semblance 30%')
axis equal, axis tight
grid
ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor = 'white';
hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
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
props = regionprops(true(size(semblance_)), semblance_, 'WeightedCentroid');
plot(props.WeightedCentroid(1)*dz, props.WeightedCentroid(2)*dx, 'gs', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'g', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);

%lgd = legend( {'true source', 'recovered source'});
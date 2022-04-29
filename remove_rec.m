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

trms_number = 2;
trms = zeros(nx, nz, trms_number);
MAPV = zeros(nx, nz, trms_number);
PAPR = zeros(nx, nz, trms_number);
for i = 2:trms_number
    filename = [num2str(sf), '_', num2str(If), '_trm_', num2str(i)];
    filenameExtention = [filename, '.mod'];
    fid = fopen(filenameExtention, 'r');
    trms(:, :, i) = fread(fid, [nx, nz], 'double');
    figure(i);
    imagesc(x, z, trms(:, :, i));
    fclose(fid);
    
        filenameMAPV = [num2str(sf), '_', num2str(If), '_MAPV_', num2str(i)];
    filenameExtentionMAPV = [filenameMAPV, '.mod'];
    fidMAPV = fopen(filenameExtentionMAPV, 'r');
    MAPV(:, :, i) = fread(fidMAPV, [nx, nz], 'double');
    fclose(fidMAPV);
    
    filenamePAPR = [num2str(sf), '_', num2str(If), '_PAPR_', num2str(i)];
    filenameExtentionPAPR = [filenamePAPR, '.mod'];
    fidPAPR = fopen(filenameExtentionPAPR, 'r');
    PAPR(:, :, i) = fread(fidPAPR, [nx, nz], 'double');
    fclose(fidPAPR);
end

avg = mean(trms(20:240, :,:), 'all');
trms(1:20,:,:)=avg;
PAPR(1:30, :, :) = 0;
semblance_ = zeros(nx, nz);
sum_max = sum(trms, 'all');
sum_max_max = max(sum(PAPR, 3),[], 'all');

for ii = 1:nx
    for jj = 1:nz
        sum_a = sum(PAPR(ii, jj, :));
        if (sum_a < 50 * sum_max_max / 100)
            semblance_(ii, jj) = 0;
        else
            num = sum(PAPR(ii, jj, :))^2;
            den = trms_number * sum(PAPR(ii, jj, :).^2); % ??
            semblance_(ii, jj) = num / den;
        end
    end
end

figure(1), clf
imagesc(x, z, semblance_);
colorbar
xlabel('x (m)')
ylabel('z (m)')
zlabel('max energy')
axis equal, axis tight

MAPV(1:30, :, :) = 0;
semblance_1 = zeros(nx, nz);
sum_max = sum(trms, 'all');
sum_max_max = max(sum(MAPV, 3),[], 'all');

for ii = 1:nx
    for jj = 1:nz
        sum_a = sum(MAPV(ii, jj, :));
        if (sum_a < 50 * sum_max_max / 100)
            semblance_1(ii, jj) = 0;
        else
            num = sum(MAPV(ii, jj, :))^2;
            den = trms_number * sum(MAPV(ii, jj, :).^2); % ??
            semblance_1(ii, jj) = num / den;
        end
    end
end

figure(3), clf
imagesc(x, z, semblance_1);
colorbar
xlabel('x (m)')
ylabel('z (m)')
zlabel('max energy')
axis equal, axis tight

figure(2);
s = surf(x, z, semblance_);
colormap('cool');
camlight;
shading interp;
set(s, 'facelighting', 'phong', 'facealpha', 0.7);
colorbar;
%axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');

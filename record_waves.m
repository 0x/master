%     (j-1/2) --- Txz-----uz-----Txz----uz-----Txz
%                  |      |      |      |      |
%                  |      |      |      |      |
%                  |      |      |      |      |
%         (j) --- ux---Txx,Tzz---ux---Txx,Tzz---ux
%                  |      |      |      |      |
%                  |      |      |      |      |
%                  |      |      |      |      |
%     (j+1/2) --- Txz-----uz-----Txz-----uz-----Txz
%                         |      |      |
%                        (i)  (i+1/2) (i+1)
% Txz(i,j) = Txz(i+1/2,j+1/2);
% ux(i,j) = ux(i+1/2,j);
% uz(i,j) = uz(i,j+1/2);

close all;
clear all;

nx = 500;
nz = 500;
dx = 4;
dz = 4;
L = nx * dx;
H = nz * dz;
receivers_number = 82;

nt_iter_stop = 4000; % Number of timesteps to compute
time_global = 1.5;
Courant = 0.9;

vp = 3000;
vs = vp / sqrt(3);
r = 2000;

sf = 300;
If = 3;
filename = [num2str(sf), '_', num2str(If), '_wavedata'];
filenameExtention = [filename, '.mod'];
fid = fopen(filenameExtention, 'w');
fprintf(fid, 'nx=%i nz=%i dx=%i dz=%i time_global=%f Courant=%f r=%f', ...
    nx, nz, dx, dz, time_global, Courant, r);
fclose(fid);

% Soure location
src_nx = floor(nx/2) + 1;
src_nz = floor(nz/2) + 1;
f0 = 30; % Hz
sc = pi * pi * f0 * f0; % Constant for source

% PML
di = 0; % Damping parameter for x
dj = 0; % Damping parameter for y
nd = 30; % Thickness of PML
pow = 2;

% create grid
x = linspace(0, L, nx);
z = linspace(0, H, nz);

% Read clatter model from file
% Compressional wave velocity [m/s]
filename = ['clatter_main_', num2str(sf), '_', num2str(If)];
filenameExtention = [filename, '.mod'];
fid = fopen(filenameExtention, 'r');
Vp = fread(fid, [nx, nz], 'double');
%Vp = ones(nx,nz)*vp;
fclose(fid);

figure(50), clf
imagesc(x, z, Vp);
colorbar
%colormap gray
xlabel('x (м)')
ylabel('z (м)')
title('\fontsize{14}Случайно-неоднородная среда')
axis equal, axis tight
Vs = Vp / sqrt(3); % Shear wave velocity [m/s]
rho = r * ones(nx, nz); % Density [kg/m^3]

miu = Vs.^2 * r;
lambda = (Vp.^2 * r - 2 * miu);

% Compute stable timestep
dt = Courant * dx / (max(max(Vp)) * sqrt(2));
nt = ceil(time_global/dt) + 1;

% Source time function
half_dur = time_global / 2.0; % Source half duration [s]
wl = min(min(Vs)) * 2.0 * half_dur; % Wavelength

% Setup initial velocity and stress profile
u = zeros(nx, nz);
ux = zeros(nx, nz);
uz = zeros(nx, nz);
v = zeros(nx, nz);
vx = zeros(nx, nz);
vz = zeros(nx, nz);

Txx = zeros(nx, nz);
Tzz = zeros(nx, nz);
Txz = zeros(nx, nz);
Txx_x = zeros(nx, nz);
Tzz_x = zeros(nx, nz);
Txz_x = zeros(nx, nz);
Txx_z = zeros(nx, nz);
Tzz_z = zeros(nx, nz);
Txz_z = zeros(nx, nz);

% Bouyancy and other parameter
b = 1 ./ rho .* ones(nx, nz);
a = dt / dx;

BU = b * a;
LAM = lambda * a;
MU = miu * a;
GAMMA = LAM + 2 * MU;

% Stations
rec_step = floor(nx/(receivers_number));
receiver = zeros(nt, receivers_number, 2); % 1 - x; 2 - z

max_energy = zeros(nx, nz, 50);
energy = zeros(nx, nz);

% plot
f = figure();
f.Position = [200 200 700 140];
tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact');
tl.Padding = 'normal';
tl.Padding = 'normal';
ax1 = nexttile(tl);
changeax1 = ones(4, 1);
changeax1(:) = true;

time = 0;
for n = 1:nt_iter_stop
    if (time >= time_global / 2 + time_global / 15 && changeax1(1))
        ax1 = nexttile(tl);
        changeax1(1) = false;
        % break;
    end
    if (time >= time_global / 2 + time_global / 8 && changeax1(2))
        ax1 = nexttile(tl);
        changeax1(2) = false;
        % break;
    end

    if (time >= time_global / 2 + time_global / 4.8 && changeax1(3))
        ax1 = nexttile(tl);
        changeax1(3) = false;
        % break;
    end
    if (time >= 1.1527)
        break;
    end
    if (time >= time_global)
        break;
    end
    %     source_term = 1 * (1 - 2 * sc * (time - half_dur)^2) * exp(-sc*(time - half_dur)^2); % Rickerra
    %     vz(src_nx, src_nz) = vz(src_nx, src_nz) + source_term * dt;
    %     ux(src_nx, src_nz) = ux(src_nx, src_nz) + source_term * dt;

    for i = 2:nx - 1
        for j = 2:nz - 1
            % PML for x
            if (i >= nx - nd)
                di = -3 * Vp(i, j) * log(0.0001) * ((+i - nx + nd) * dx)^pow / (2 * (nd * dx)^3);
            else
                di = 0;
            end
            % PML for z
            if (j >= nz - nd)
                dj = -3 * Vp(i, j) * log(0.0001) * ((+j - nz + nd) * dz)^pow / (2 * (nd * dz)^3);
            elseif (j <= nd)
                dj = -3 * Vp(i, j) * log(0.0001) * (nd * dz - j * dz)^pow / (2 * (nd * dz)^3);
            else
                dj = 0;
            end

            ux(i, j) = (1 - dt * di / 2) * ux(i, j) / (1 + dt * di / 2) + BU(i, j) * (Txx(i + 1, j) - Txx(i, j)) / (1 + dt * di / 2);
            uz(i, j) = (1 - dt * dj / 2) * uz(i, j) / (1 + dt * dj / 2) + BU(i, j) * (Txz(i, j) - Txz(i, j - 1)) / (1 + dt * dj / 2);
            u(i, j) = ux(i, j) + uz(i, j);

            vx(i, j) = (1 - dt * di / 2) * vx(i, j) / (1 + dt * di / 2) + BU(i, j) * (Txz(i, j) - Txz(i - 1,j)) / (1 + dt * di / 2);
            vz(i, j) = (1 - dt * dj / 2) * vz(i, j) / (1 + dt * dj / 2) + BU(i, j) * (Tzz(i, j + 1) - Tzz(i, j)) / (1 + dt * dj / 2);
            v(i, j) = vx(i, j) + vz(i, j);

        end
    end

    source_term = 1 * (1 - 2 * sc * (time - half_dur)^2) * exp(-sc*(time - half_dur)^2); % Rickerra
    Txx_x(src_nx, src_nz) = Txx_x(src_nx, src_nz) + source_term * dt;
    Tzz_z(src_nx, src_nz) = Tzz_z(src_nx, src_nz) + source_term * dt;

    for i = 2:nx - 1
        for j = 2:nz - 1
            % PML for x
            if (i >= nx - nd)
                di = -3 * Vp(i, j) * log(0.0001) * ((+i - nx + nd) * dx)^pow / (2 * (nd * dx)^3);
            else
                di = 0;
            end
            % PML for z
            if (j >= nz - nd)
                dj = -3 * Vp(i, j) * log(0.0001) * ((+j - nz + nd) * dz)^pow / (2 * (nd * dz)^3);
            elseif (j <= nd)
                dj = -3 * Vp(i, j) * log(0.0001) * (nd * dz - j * dz)^pow / (2 * (nd * dz)^3);
            else
                dj = 0;
            end

            Txx_x(i, j) = (1 - dt * di / 2) * Txx_x(i, j) / (1 + dt * di / 2) + GAMMA(i, j) * (u(i, j) - u(i - 1, j)) / (1 + dt * di / 2);
            Txx_z(i, j) = (1 - dt * dj / 2) * Txx_z(i, j) / (1 + dt * dj / 2) + LAM(i, j) * (v(i, j) - v(i, j - 1)) / (1 + dt * dj / 2);
            Txx(i, j) = Txx_x(i, j) + Txx_z(i, j);

            Tzz_x(i, j) = (1 - dt * di / 2) * Tzz_x(i, j) / (1 + dt * di / 2) + LAM(i, j) * (u(i, j) - u(i - 1, j)) / (1 + dt * di / 2);
            Tzz_z(i, j) = (1 - dt * dj / 2) * Tzz_z(i, j) / (1 + dt * dj / 2) + GAMMA(i, j) * (v(i, j) - v(i, j - 1)) / (1 + dt * dj / 2);
            Tzz(i, j) = Tzz_x(i, j) + Tzz_z(i, j);

            Txz_x(i, j) = (1 - dt * di / 2) * Txz_x(i, j) / (1 + dt * di / 2) + MU(i, j) * (v(i + 1, j) - v(i, j)) / (1 + dt * di / 2);
            Txz_z(i, j) = (1 - dt * dj / 2) * Txz_z(i, j) / (1 + dt * dj / 2) + MU(i, j) * (u(i, j + 1) - u(i, j)) / (1 + dt * dj / 2);
            Txz(i, j) = Txz_x(i, j) + Txz_z(i, j);
        end
    end

    rec_coord = 0;
    for i = 1:receivers_number
        rec_coord = rec_coord + rec_step;
        receiver(n, i, 1) = Txx(2, rec_coord);
        receiver(n, i, 2) = Tzz(2, rec_coord);
    end

    time = time + dt;

    % Plot solution every 50 timesteps
    if (mod(n, 50) == 0)
        imagesc(ax1, x, z, sqrt(Txx.^2 + Tzz.^2));
        colorbar;
        axis equal;
        axis tight;
        %colormap gray;
        xlabel(tl, 'x (м)')
        ylabel(tl, 'z (м)')
        %zlabel(tl,'Pressure (Pa)')
        title(['\fontsize{14}', num2str(time), ' сек.'])
        axis equal, axis tight
hold on
                plot(src_nz*dz, src_nx*dx, '.', ...
                    'MarkerSize', 10, ...
                    'MarkerEdgeColor', 'm', ...
                    'MarkerFaceColor', [0.5, 0.5, 0.5]);
                drawnow
    end
    


end
legend('\fontsize{12}источник')

figure();
imagesc(x, z, sqrt(Txx.^2 + Tzz.^2));
colorbar;
%colormap gray;
xlabel('x (м)')
ylabel('z (м)')
title(['\fontsize{14}', num2str(time), ' сек.'])
axis equal, axis tight

hold on
plot(src_nz*dz, src_nx*dx, 'ms', ...
    'LineWidth', 2, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'm', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
hold on

x_recv = zeros(receivers_number, 1);
y_recv = zeros(receivers_number, 1);

rec_coord = 0;
for i = 1:receivers_number
    rec_coord = rec_coord + rec_step;
    x_recv(i) = rec_coord * dz;
    y_recv(i) = 2 * dx;

end
plot(x_recv, y_recv, 'cv', ...
    'LineWidth', 2, ...
    'MarkerSize', 7, ...
    'MarkerEdgeColor', 'c', ...
    'MarkerFaceColor', [0.5, 0.5, 0.5])
legend('\fontsize{12}источник', '\fontsize{12}приемник')

% for i = 1:receivers_number
%     filename = [num2str(sf), '_', num2str(If), '_receiver82_x_', num2str(i)];
%     filenameExtention = [filename, '.mod'];
%     fid = fopen(filenameExtention, 'w');
%     fwrite(fid, receiver(:, i, 1), 'double');
%     fclose(fid);
%     filename = [num2str(sf), '_', num2str(If), '_receiver82_z_', num2str(i)];
%     filenameExtention = [filename, '.mod'];
%     fid = fopen(filenameExtention, 'w');
%     fwrite(fid, receiver(:, i, 2), 'double');
%     fclose(fid);
% end

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

nt_iter_stop = 3000; % Number of timesteps to compute
time_global = 1.5;

vp = 3000;
vs = vp / sqrt(3);
r = 2000;

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
sf = 150;
If = 2;
filename = ['clatter_main_', num2str(sf), '_', num2str(If)];
filenameExtention = [filename, '.mod'];
fid = fopen(filenameExtention, 'r');
Vp = fread(fid, [nx, nz], 'double');
%Vp=ones(nx,nz)*vp;
fclose(fid);

figure(50), clf
imagesc(x, z, Vp);
colorbar
%colormap gray
xlabel('x (m)')
ylabel('z (m)')
title('clutter')
axis equal, axis tight

Vs = Vp / sqrt(3); % Shear wave velocity [m/s]
rho = r * ones(nx, nz); % Density [kg/m^3]

miu = Vs.^2 * r;
lambda = (Vp.^2 * r - 2 * miu);

% Compute stable timestep
dt = 0.9 * dx / (max(max(Vp)) * sqrt(2));

nt = ceil(time_global/dt) + 1;

% Source time function
half_dur = time_global / 2.0; % Source half duration [s]
nt_half = ceil(nt/2);
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
rec_step = floor(nx/(receivers_number + 1));
receiver = zeros(nt, receivers_number, 2); % 1 - x; 2 - z

for i = 1:receivers_number
    filename = [num2str(sf), '_', num2str(If), '_receiver82_x_', num2str(i)];
    filenameExtention = [filename, '.mod'];
    fid = fopen(filenameExtention, 'r');
    receiver(:, i, 1) = fread(fid, nt, 'double');
    fclose(fid);

    filename = [num2str(sf), '_', num2str(If), '_receiver82_z_', num2str(i)];
    filenameExtention = [filename, '.mod'];
    fid = fopen(filenameExtention, 'r');
    receiver(:, i, 2) = fread(fid, nt, 'double');
    fclose(fid);
end

impls_number = 1;
imax = zeros(impls_number);
jmax = zeros(impls_number);

max_energy = zeros(nx, nz, impls_number);
energy = zeros(nx, nz, impls_number);
MAPV = zeros(nx, nz, impls_number);
PAPR = zeros(nx, nz, impls_number);
P_total = zeros(nx, nz);
P_max = zeros(nx, nz);

% plot
f = figure();
f.Position = [200 200 700 140];
tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact');
tl.Padding = 'normal';
ax1 = nexttile(tl);
changeax1 = ones(4, 1);
changeax1(:) = true;
disp('1ok');

for impl_number = 1:impls_number
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

    time = 0;

    %filename = ['clt_generated' , num2str(layer), '_', num2str(impl_number), '_', num2str(sf), '_', num2str(If)];
    filename = ['clatter_main_', num2str(sf), '_', num2str(If)];
    filenameExtention = [filename, '.mod'];
    fid = fopen(filenameExtention, 'r');
    %Vp = fread(fid, [nx, nz], 'double');
    Vp=ones(nx,nz)*vp;
    fclose(fid);

    % recalculate all params
    Vs = Vp / sqrt(3); % Shear wave velocity [m/s]
    rho = r * ones(nx, nz); % Density [kg/m^3]

    miu = Vs.^2 * r;
    lambda = (Vp.^2 * r - 2 * miu);

    % Bouyancy and other parameter
    b = 1 ./ rho .* ones(nx, nz);
    a = dt / dx;

    BU = b * a;
    LAM = lambda * a;
    MU = miu * a;
    GAMMA = LAM + 2 * MU;

    for n = nt:-1:nt_half
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

        if (n <= nt)
            rec_coord = 0;
            for i = 1:receivers_number
                rec_coord = rec_coord + rec_step;
                Txx_x(2, rec_coord) = Txx_x(2, rec_coord) + receiver(n, i, 1);
                Tzz_z(2, rec_coord) = Tzz_z(2, rec_coord) + receiver(n, i, 2);
            end
        end

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

        % Calculate energy

        nu = 0.5 * lambda ./ (lambda + miu);
        E = 2.0 * miu .* (1.0 + nu);
        % Compute total field from split components
        epsilon_xx = (Txx - nu .* Tzz) ./ E;
        epsilon_zz = (Tzz - nu .* Txx) ./ E;
        epsilon_xz = 2.0 * (1.0 + nu) .* Txz ./ E;

        energy(:, :, impl_number) = 0.50 * (0.5 * (r + r) .* u.^2 ...
            +0.5 * (r + r) .* v.^2) + ...
            0.50 * (epsilon_xx .* Txx + epsilon_zz .* Tzz + 2 * epsilon_xz .* Txz);

        max_energy(:, :, impl_number) = (max_energy(:, :, impl_number) + energy(:, :, impl_number));

        MAPV(:, :, impl_number) = max(MAPV(:, :, impl_number), abs(energy(:, :, impl_number)));
        P_max = max(P_max, energy(:, :, impl_number));
        P_total = P_total + energy(:, :, impl_number);

        if (n == nt - 620 && changeax1(1))
            ax1 = nexttile(tl);
            changeax1(1) = false;
        end
        if (n == nt - 769 && changeax1(2))
            ax1 = nexttile(tl);
            changeax1(2) = false;
        end

        if (n == nt - 918 && changeax1(3))
            ax1 = nexttile(tl);
            changeax1(3) = false;
        end
        if (n == nt - 1068)
            %break;
        end


        % Plot solution every 50 timesteps
        if (true)
            if (n == nt - 917 || n == nt - 767 || n == nt - 618 || n == nt - 1066)
                %figure(661), clf
                imagesc(ax1, x, z, sqrt(Txx.^2 + Tzz.^2));
                colorbar
                xlabel(tl, 'x (m)')
                ylabel(tl, 'z (m)')
                title(['\fontsize{14}', num2str(time), ' sec.'])
                axis equal, axis tight
                hold on
                plot(src_nz*dz, src_nx*dx, '.', ...
                    'MarkerSize', 10, ...
                    'MarkerEdgeColor', 'm', ...
                    'MarkerFaceColor', [0.5, 0.5, 0.5]);
                legend('source');
                drawnow
                %             figure(991), clf
                %             imagesc(x,z,energy(:, :, impl_number));
                %             colorbar
                %             xlabel('x (m)')
                %             ylabel('z (m)')
                %             zlabel('Pressure (Pa)')
                %             title(['Time = ', num2str(time), ' sec'])
                %             axis equal, axis tight
                %             drawnow
            end
        end
        %legend('\fontsize{12}источник')


        time = time + dt;
    end

    PAPR(:, :, impl_number) = (P_max.^2) ./ (mean(energy(:, :, impl_number), 'all'));
    if (true)
        disp('ok');
        %filenamePAPR = [num2str(sf), '_', num2str(If), '_PAPR_generated' , num2str(layer), '_', num2str(impl_number)];
        filenamePAPR = [num2str(sf), '_', num2str(If), '_PAPR_het_', num2str(impl_number)];
        filenameExtentionPAPR = [filenamePAPR, '.mod'];
        fidPAPR = fopen(filenameExtentionPAPR, 'w');
        fwrite(fidPAPR, PAPR(:, :, impl_number), 'double');
        fclose(fidPAPR);

        %filenameLastEnergy = [num2str(sf), '_', num2str(If), '_energy_generated' , num2str(layer), '_', num2str(impl_number)];
        filenameLastEnergy = [num2str(sf), '_', num2str(If), '_energy_het_', num2str(impl_number)];
        filenameExtentionLastEnergy = [filenameLastEnergy, '.mod'];
        fidLastEnergy = fopen(filenameExtentionLastEnergy, 'w');
        fwrite(fidLastEnergy, energy(:, :, impl_number), 'double');
        fclose(fidLastEnergy);

                    filename = [num2str(sf), '_', num2str(If), '_trm_het_', num2str(impl_number)];
                    filenameExtention = [filename, '.mod'];
                    fid = fopen(filenameExtention, 'w');
                    fwrite(fid, max_energy(:, :, impl_number), 'double');
                    fclose(fid);


                filename = [num2str(sf), '_', num2str(If), '_lastwaveTxx_het_', num2str(impl_number)];
                filenameExtention = [filename, '.mod'];
                fid = fopen(filenameExtention, 'w');
                fwrite(fid, Txx, 'double');
                fclose(fid);
        
                filename = [num2str(sf), '_', num2str(If), '_lastwaveTzz_het_', num2str(impl_number)];
                filenameExtention = [filename, '.mod'];
                fidTzz = fopen(filenameExtention, 'w');
                fwrite(fidTzz, Tzz, 'double');
                fclose(fidTzz);
    end
    max_value = max_energy(1, 1, impl_number);
    for i = 30:nx
        for j = nd:nz
            if (max_value < max_energy(i, j, impl_number))
                max_value = max_energy(i, j, impl_number);
                imax(impl_number) = i;
                jmax(impl_number) = j;
            end
        end
    end

%     figure(impl_number), clf
%     imagesc(x, z, max_energy(:, :, impl_number));
%     colorbar
%     xlabel('x (м)')
%     ylabel('z (м)')
%     zlabel('max energy')
%     title(['\fontsize{14}', num2str(time), ' сек.'])
%     axis equal, axis tight
%     hold on;
%     plot(jmax(impl_number)*dz, imax(impl_number)*dx, 'gs', ...
%         'LineWidth', 2, ...
%         'MarkerSize', 10, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', [0.5, 0.5, 0.5]);
% 
%     figure();
%     s = surf(x, z, max_energy(:, :, impl_number));
%     colormap('cool');
%     camlight;
%     shading interp;
%     set(s, 'facelighting', 'phong', 'facealpha', 0.7);
%     colorbar;
%     %axis equal;
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
end

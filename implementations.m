clear all;
close all;

implementations_number = 200;
Nx = 500; % grid points
Nz = 500;

hx = 4; % grid size
hy = 4;

% Matematicheskoe ozhidanie
mf = 3000; % middle velocity

% Lx=Nx*hx;     % meters
% Lz=Nz*hz;     % meters

% Don't change this
Ny = 5;
hz = 4;

% Standartnoe otklonenie
sf = 300;
% Dlina korrelyacii
If = 25;

lf = 2 * If / sqrt(pi);

Ng = 1000;
C = sf / sqrt(Ng);

for iter = 12:implementations_number
    Yf = zeros(Nx, Ny, Nz);
    for ig = 1:Ng

        fi = 2.0 * pi * rand;

        % Gaussian spectrum
        flag = true;
        while (flag)
            k = 100 * rand;
            d = k * k * exp(-0.5*k*k);
            if (rand * 2 * exp(-1) < d)
                flag = false;
            end
        end

        k = sqrt(2) * k / lf;

        theta = acos(1-2*rand);

        V(1) = k * sin(fi) * sin(theta);
        V(2) = k * cos(fi) * sin(theta);
        V(3) = k * cos(theta);

        a = randn;
        b = randn;

        skx = zeros(Nx, Ny, Nz);

        for ix = 1:Nx
            for iy = 1:Ny
                for iz = 1:Nz
                    skx(ix, iy, iz) = hx * (ix - 0.5) * V(1) + hy * (iy - 0.5) * V(2) + hz * (iz - 0.5) * V(3);
                end
            end
        end
        Yf = Yf + a * cos(skx) + b * sin(skx);

    end
    Yf = C * Yf;

    Fxy = zeros(Nx, Ny);
    Fxz = zeros(Nx, Nz);
    Fyz = zeros(Ny, Nz);

    Fxy(:, :) = Yf(:, :, round(Nz / 2));
    Fxz(:, :) = Yf(:, round(Ny / 2), :);
    Fyz(:, :) = Yf(round(Nx / 2), :, :);

    filename = ['clt_', num2str(iter), '_', num2str(sf), '_', num2str(If)];
    disp(filename)

    f1 = fopen([filename, '.mod'], 'w');
    fwrite(f1, mf+Fxz, 'double');
    fclose(f1);
    figure();
    imagesc((mf+Fxz)');close
    title(['Standard deviation=', num2str(sf), '; correlation length=', num2str(If)])
    colorbar;
    axis equal, axis tight
end

% figure;
% imagesc((mf+Fxz)');
% title(['Standard deviation=', num2str(sf), '; correlation length=', num2str(If)])
% colorbar;
% axis equal, axis tight

% f1 = fopen([filename, '.modinfo'], 'w');
% fprintf(f1, 'n1=%i n2=%i d1=%i d2=%i in=./%s esize=4 data_format="native_float"', ...
%     Nz, Nx, hz, hx, [filename, '.mod']);
% fclose(f1);

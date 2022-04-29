

%Name1 ='Model_clatter_Vp_S_20_I_4.mod';
Name1 ='clatter_main_300_3.mod';
%Name1 ='wave.mod';
fid = fopen( Name1 ,'r');

Nx=500; % grid points
Nz=500;

hx=4; % grid size
hz=4;

Lx=(Nx-1)*hx;     % meters
Lz=(Nz-1)*hz;     % meters
% create grid
x = linspace(0, L, nx);
z = linspace(0, H, nz);
f = figure();
f.Position = [200 200 380 280];
u=fread(fid,[Nx,Nz],'float64');
imagesc(x,z,u);
%grid
xlabel('x (m)')
ylabel('z (m)')
title(['\fontsize{14}Randomly inhomogeneous medium']);
axis equal;
axis tight;
colorbar;



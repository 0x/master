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

nx = 201;
nz = 201;

dx = 20;
dz = 20;

L = nx*dx;
H = nz*dz;

nt_iter_stop = 3000;  % Number of timesteps to compute
time_global = 1.3;

vp = 5500;
vs = vp/sqrt(3);
r = 2600;

% Source location
src_nx = floor(nz/4)+1;
src_nz = floor(nz/2)+1;
f0=15;  % Hz
ac=pi*pi*f0*f0;  % Constant for source

% PML
di=0;  % Damping parameter for x
dj=0;  % Damping parameter for y
nd=20;  % Thickness of PML
pow = 2;

Vp = vp*ones(nx,nz);  % Compressional wave velocity [m/s]
Vs = vs*ones(nx,nz);  % Shear wave velocity [m/s]
rho = r*ones(nx,nz);  % Density [kg/m^3]

miu = Vs.^2*r;
lambda = (Vp.^2*r-2*miu);

% Compute stable timestep
dt = 0.9*dx/(vp*sqrt(2));
nt = ceil(time_global/dt)+1;

% Source time function
half_dur = 0.2;  % Source half duration [s]
wl = min(min(Vs))*2*half_dur;  % Wavelength

% Setup initial velocity and stress profile
u = zeros(nx,nz);
ux = zeros(nx,nz);
uz = zeros(nx,nz);
v = zeros(nx,nz);
vx = zeros(nx,nz);
vz = zeros(nx,nz);

Txx = zeros(nx,nz);
Tzz = zeros(nx,nz);
Txz = zeros(nx,nz);
Txx_x = zeros(nx,nz);
Tzz_x = zeros(nx,nz);
Txz_x = zeros(nx,nz);
Txx_z = zeros(nx,nz);
Tzz_z = zeros(nx,nz);
Txz_z = zeros(nx,nz);

% Bouyancy and other parameter
b = 1./rho.*ones(nx,nz);
a = dt/dx;

BU = b*a;
LAM = lambda*a;
MU = miu*a;
GAMMA = LAM + 2*MU;

% create grid
x=linspace(0,L,nx);
z=linspace(0,H,nz);

% Stations
receiver1 = zeros(nt,1);
receiver2 = zeros(nt,1);
receiver3 = zeros(nt,1);
receiver4 = zeros(nt,1);
receiver5 = zeros(nt,1);
receiver6 = zeros(nt,1);
receiver7 = zeros(nt,1);
receiver8 = zeros(nt,1);
receiver9 = zeros(nt,1);
max_energy=zeros(nx, nz);
energy=zeros(nx, nz);

time = 0;
for n=1:nt_iter_stop
    if (time >= time_global)
        break;
    end
    
    for i=2:nx-1
        for j=2:nz-1
            % PML for x
            if (i>=nx-nd)
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)^pow/(2*(nd*dx)^3);
            else
                di=0;
            end
            % PML for z
            if (j>=nz-nd)
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)^pow/(2*(nd*dz)^3);
            elseif (j<=nd)
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)^pow/(2*(nd*dz)^3);
            else
                dj=0;
            end
            
            ux(i,j) = (1-dt*di/2)*ux(i,j)/(1+dt*di/2) + BU(i,j) * (Txx(i+1,j)-Txx(i,j))/(1+dt*di/2);
            uz(i,j) = (1-dt*dj/2)*uz(i,j)/(1+dt*dj/2) + BU(i,j) * (Txz(i,j)-Txz(i,j-1))/(1+dt*dj/2);
            u(i,j) = ux(i,j) + uz(i,j);
            
            vx(i,j) = (1-dt*di/2)*vx(i,j)/(1+dt*di/2) + BU(i,j) * (Txz(i,j) - Txz(i - 1,j))/(1+dt*di/2);
            vz(i,j) = (1-dt*dj/2)*vz(i,j)/(1+dt*dj/2) + BU(i,j) * (Tzz(i,j + 1) -  Tzz(i,j))/(1+dt*dj/2);
            v(i,j) = vx(i,j) + vz(i,j);
            
            if (i == 2 && j == 40)
                receiver1(n) = u(i,j);
            elseif (i == 2 && j == 60)
                receiver2(n) = u(i,j);
            elseif (i == 2 && j == 80)
                receiver3(n) = u(i,j);
            elseif (i == 2 && j == 100)
                receiver4(n) = u(i,j);
            elseif (i == 2 && j == 120)
                receiver5(n) = u(i,j);
            elseif (i == 2 && j == 140)
                receiver6(n) = u(i,j);
            elseif (i == 2 && j == 160)
                receiver7(n) = u(i,j);
            elseif (i == 2 && j == 180)
                receiver8(n) = u(i,j);
            elseif (i == 2 && j == 20)
                receiver9(n) = u(i,j);
            end
            
        end
    end
       
    source_term=1*(1-2*ac*(time-half_dur)^2)*exp(-ac*(time-half_dur)^2); % Rickerra
    Tzz_z(src_nx,src_nz) =  Tzz_z(src_nx,src_nz)+source_term*dt;
    %Txx_x(src_nx,src_nz) =  Txx_x(src_nx,src_nz)+source_term*dt;
    
    for i=2:nx-1
        for j=2:nz-1
            % PML for x
            if (i>=nx-nd)
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)^pow/(2*(nd*dx)^3);
            else
                di=0;
            end
            % PML for z
            if (j>=nz-nd)
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)^pow/(2*(nd*dz)^3);
            elseif (j<=nd)
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)^pow/(2*(nd*dz)^3);
            else
                dj=0;
            end
            
            Txx_x(i,j) = (1-dt*di/2)*Txx_x(i,j)/(1+dt*di/2) + GAMMA(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2);
            Txx_z(i,j) = (1-dt*dj/2)*Txx_z(i,j)/(1+dt*dj/2) + LAM(i,j)*(v(i,j)-v(i,j-1))/(1+dt*dj/2);
            Txx(i,j) = Txx_x(i,j) + Txx_z(i,j);
            
            Tzz_x(i,j) = (1-dt*di/2)*Tzz_x(i,j)/(1+dt*di/2) + LAM(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2);
            Tzz_z(i,j) = (1-dt*dj/2)*Tzz_z(i,j)/(1+dt*dj/2) + GAMMA(i,j)*(v(i,j)- v(i,j-1))/(1+dt*dj/2);
            Tzz(i,j) = Tzz_x(i,j) + Tzz_z(i,j);
            
            Txz_x(i,j) = (1-dt*di/2)*Txz_x(i,j)/(1+dt*di/2) + MU(i,j) * (v(i+1,j)- v(i,j))/(1+dt*di/2);
            Txz_z(i,j) = (1-dt*dj/2)*Txz_z(i,j)/(1+dt*dj/2) + MU(i,j)*(u(i,j+1)-u(i,j))/(1+dt*dj/2);
            Txz(i,j) = Txz_x(i,j) + Txz_z(i,j);
        end
    end
    
    time = time+dt;
    
    %     % Plot solution every 50 timesteps
    %     if (mod(n,5)==0)
    %         figure(1), clf
    %         imagesc(x,z,uz);
    %         colorbar
    %         xlabel('x [m]')
    %         ylabel('z [m]')
    %         zlabel('Pressure [Pa]')
    %         title(['Time = ',num2str(time),' sec'])
    %         axis equal, axis tight
    %         drawnow
    %     end
end

figure(1), clf
imagesc(x,z,v);
colorbar
xlabel('x [m]')
ylabel('z [m]')
zlabel('Pressure [Pa]')
title(['Time = ',num2str(time),' sec'])
axis equal, axis tight
hold on
plot(src_nz*dz, src_nx*dx,'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5])

u = zeros(nx,nz);
ux = zeros(nx,nz);
uz = zeros(nx,nz);
v = zeros(nx,nz);
vx = zeros(nx,nz);
vz = zeros(nx,nz);

Txx = zeros(nx,nz);
Tzz = zeros(nx,nz);
Txz = zeros(nx,nz);
Txx_x = zeros(nx,nz);
Tzz_x = zeros(nx,nz);
Txz_x = zeros(nx,nz);
Txx_z = zeros(nx,nz);
Tzz_z = zeros(nx,nz);
Txz_z = zeros(nx,nz);

time = 0;
for n=nt:-1:1
    for i=2:nx-1
        for j=2:nz-1
            % PML for x
            if (i>=nx-nd)
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)^pow/(2*(nd*dx)^3);
            else
                di=0;
            end
            % PML for z
            if (j>=nz-nd)
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)^pow/(2*(nd*dz)^3);
            elseif (j<=nd)
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)^pow/(2*(nd*dz)^3);
            else
                dj=0;
            end
            
            ux(i,j) = (1-dt*di/2)*ux(i,j)/(1+dt*di/2) + BU(i,j) * (Txx(i+1,j)-Txx(i,j))/(1+dt*di/2);
            uz(i,j) = (1-dt*dj/2)*uz(i,j)/(1+dt*dj/2) + BU(i,j) * (Txz(i,j)-Txz(i,j-1))/(1+dt*dj/2);
            u(i,j) = ux(i,j) + uz(i,j);
            
            vx(i,j) = (1-dt*di/2)*vx(i,j)/(1+dt*di/2) + BU(i,j) * (Txz(i,j) - Txz(i - 1,j))/(1+dt*di/2);
            vz(i,j) = (1-dt*dj/2)*vz(i,j)/(1+dt*dj/2) + BU(i,j) * (Tzz(i,j + 1) -  Tzz(i,j))/(1+dt*dj/2);
            v(i,j) = vx(i,j) + vz(i,j); 
        end
    end
    
    if (n <= nt)
        Tzz_z(2,40)=Tzz_z(2,40)+receiver1(n);
        Tzz_z(2,60)=Tzz_z(2,60)+receiver2(n);
        Tzz_z(2,80)=Tzz_z(2,80)+receiver3(n);
        Tzz_z(2,100)=Tzz_z(2,100)+receiver4(n);
        Tzz_z(2,120)=Tzz_z(2,120)+receiver5(n);
        Tzz_z(2,140)=Tzz_z(2,140)+receiver6(n);
        Tzz_z(2,160)=Tzz_z(2,160)+receiver7(n);
        Tzz_z(2,180)=Tzz_z(2,180)+receiver8(n);
        Tzz_z(2,20)=Tzz_z(2,20)+receiver9(n);
    end
    
    for i=2:nx-1
        for j=2:nz-1
            % PML for x
            if (i>=nx-nd)
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)^pow/(2*(nd*dx)^3);
            else
                di=0;
            end
            % PML for z
            if (j>=nz-nd)
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)^pow/(2*(nd*dz)^3);
            elseif (j<=nd)
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)^pow/(2*(nd*dz)^3);
            else
                dj=0;
            end
            
            Txx_x(i,j) = (1-dt*di/2)*Txx_x(i,j)/(1+dt*di/2) + GAMMA(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2);
            Txx_z(i,j) = (1-dt*dj/2)*Txx_z(i,j)/(1+dt*dj/2) + LAM(i,j)*(v(i,j)-v(i,j-1))/(1+dt*dj/2);
            Txx(i,j) = Txx_x(i,j) + Txx_z(i,j);
            
            Tzz_x(i,j) = (1-dt*di/2)*Tzz_x(i,j)/(1+dt*di/2) + LAM(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2);
            Tzz_z(i,j) = (1-dt*dj/2)*Tzz_z(i,j)/(1+dt*dj/2) + GAMMA(i,j)*(v(i,j)- v(i,j-1))/(1+dt*dj/2);
            Tzz(i,j) = Tzz_x(i,j) + Tzz_z(i,j);
            
            Txz_x(i,j) = (1-dt*di/2)*Txz_x(i,j)/(1+dt*di/2) + MU(i,j) * (v(i+1,j)- v(i,j))/(1+dt*di/2);
            Txz_z(i,j) = (1-dt*dj/2)*Txz_z(i,j)/(1+dt*dj/2) + MU(i,j)*(u(i,j+1)-u(i,j))/(1+dt*dj/2);
            Txz(i,j) = Txz_x(i,j) + Txz_z(i,j);
        end
    end
    
    % Calculate energy
    for ii=1:nx
        for jj=1:nz
            nu = 0.5 * lambda(ii,jj)/(lambda(ii,jj) + miu(ii,jj));
            E = 2.0*miu(ii,jj)*(1.0+ nu);
            % Compute total field from split components
            epsilon_xx = ( Txx(ii,jj) - nu*Tzz(ii,jj) )/E;
            epsilon_zz = ( Tzz(ii,jj) - nu*Txx(ii,jj) )/E;
            epsilon_xz = 2.0*(1.0+ nu)*Txz(ii,jj) / E;
            
            energy(ii, jj) = 0.50 * ( 0.5*(r+r)* u(ii,jj)^2 ...
                + 0.5*(r+r)*v(ii,jj)^2) + ...
                0.50 * ( epsilon_xx * Txx(ii,jj) + epsilon_zz * Tzz(ii,jj) + epsilon_xz * Txz(ii,jj));
            
            max_energy(ii, jj) = (max_energy(ii, jj)+ energy(ii, jj));
        end
    end

    % Plot solution every 50 timesteps
    if (mod(n,5)==0)
        figure(2), clf
        imagesc(x,z,uz);
        colorbar
        xlabel('x [m]')
        ylabel('z [m]')
        zlabel('Pressure [Pa]')
        title(['Time = ',num2str(time),' sec'])
        axis equal, axis tight
        drawnow
    end
    
    time = time+dt;
end

max_value = max_energy(1,1);
imax = 0;
jmax = 0;
for i=10:nx
    for j=nd:nz
        if (max_value < max_energy(i,j))
            max_value = max_energy(i,j);
            imax = i;
            jmax = j;
        end
    end
end

figure(3), clf
imagesc(x,z,max_energy);
colorbar
xlabel('x [m]')
ylabel('z [m]')
zlabel('maxEnergy')
title(['Time = ',num2str(time),' sec'])
axis equal, axis tight
hold on;
plot(jmax*dz, imax*dx,'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5])


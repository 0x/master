clear; close('all')

nx=200;  % Number points of x
nz=200;  % Number points of z

dx=20;  % 20 m
dz=dx;
T=1.1;  % Total time

isrc = 50;
jsrc = 100;

vp=zeros(nx,nz);
vs=zeros(nx,nz);

% Read vp clutter model from file
filename = 'Model_clatter_Vp_S_20_I_4';
filenameExtention =[filename '.mod'];
fid = fopen(filenameExtention ,'r');
vp=fread(fid,[nx,nz],'double');
imagesc(vp);
title('clutter; Nx=Nz=200; hx=hz=20')
colorbar;

ro=2000;
%vp(:)=5500;
%vs=vp/sqrt(3);
vs=0;

lambda=ro*(vp.^2-2*vs.^2);  % 2.621666666666665e+10
nyu=ro*vs.^2;  % 2.621666666666665e+10

CFL=0.5;  % CFL=vp*dt/dx
dt=CFL*min(dx,dz)/(max(max(vp/2))*sqrt(2.0));

t0 = T/2;
f0=15;  % Hz
a=pi*pi*f0*f0;  % Constant for source

Lx=nx*dx;  % Size of domain
Lz=nz*dz;

x=linspace(0,Lx,nx);
z=linspace(0,Lz,nz);

% V - Vz
% U - Vx

% V
v=zeros(nx,nz);     % v(n), optional
vx=zeros(nx,nz);    % vx(n)
vz=zeros(nx,nz);    % vz(n)
vxm1=zeros(nx,nz);  % vx(n-1)
vzm1=zeros(nx,nz);  % vz(n-1)
% U
u=zeros(nx,nz);     % u(n), optional
ux=zeros(nx,nz);    % ux(n)
uz=zeros(nx,nz);    % uz(n)
uxm1=zeros(nx,nz);  % ux(n-1)
uzm1=zeros(nx,nz);  % uz(n-1)
% tau(xx)
tau1=zeros(nx,nz);     % tauxx(n), optional
tau1x=zeros(nx,nz);    % tauxx(x)(n)
tau1z=zeros(nx,nz);    % tauxx(z)(n)
tau1xp1=zeros(nx,nz);  % tauxx(x)(n+1)
tau1zp1=zeros(nx,nz);  % tauxx(z)(n+1)
% tau(xz)
tau2=zeros(nx,nz);     % tauxz(n), optional
tau2x=zeros(nx,nz);    % tauxz(x)(n)
tau2z=zeros(nx,nz);    % tauxz(z)(n)
tau2xp1=zeros(nx,nz);  % tauxz(x)(n+1)
tau2zp1=zeros(nx,nz);  % tauxz(z)(n+1)
% tau(zz)
tau3=zeros(nx,nz);     % tauzz(n), optional
tau3x=zeros(nx,nz);    % tauzz(x)(n)
tau3z=zeros(nx,nz);    % tauzz(z)(n)
tau3xp1=zeros(nx,nz);  % tauzz(x)(n+1)
tau3zp1=zeros(nx,nz);  % tauzz(z)(n+1)

B=zeros(nx,nz);  % Buoyancy
L=zeros(nx,nz);  % Lame coef
M=zeros(nx,nz);  % Lame coef

% PML
di=0;   % Damping parameter for x
dj=0;   % Damping parameter for y
nd=5;  % Thickness of PML
pow = 2;

% Make constant for now
B(:)=1/ro;
L(:)=lambda;
M(:)=nyu;
iter =0;
% Scheme
t=0;

nt = ceil(T/dt)+1;
receiver1 = zeros(nt,1);
receiver2 = zeros(nt,1);
receiver3 = zeros(nt,1);
receiver4 = zeros(nt,1);
receiver5 = zeros(nt,1);
receiver6 = zeros(nt,1);
receiver7 = zeros(nt,1);
receiver8 = zeros(nt,1);
receiver9 = zeros(nt,1);

iter_t = 1;

while (t<T)
    t=t+dt;
    
    % Source
    source_term=1*(1-2*a*(t-T/2)^2)*exp(-a*(t-T/2)^2); % Rickerra
    tau3(isrc,jsrc)=tau3(isrc,jsrc)+source_term*dt;
    %tau1(101,101)=tau1(101,101)+source_term*dt;
    %disp(source_term);
    
    for i=2:nx-1
        for j=2:nz-1
            if (i == 2 && j == 40)
                receiver1(iter_t) = u(i,j);
            elseif (i == 2 && j == 60)
                receiver2(iter_t) = u(i,j);
            elseif (i == 2 && j == 80)
                receiver3(iter_t) = u(i,j);
            elseif (i == 2 && j == 100)
                receiver4(iter_t) = u(i,j);
            elseif (i == 2 && j == 120)
                receiver5(iter_t) = u(i,j);
            elseif (i == 2 && j == 140)
                receiver6(iter_t) = u(i,j);
            elseif (i == 2 && j == 160)
                receiver7(iter_t) = u(i,j);
            elseif (i == 2 && j == 180)
                receiver8(iter_t) = u(i,j);
             elseif (i == 2 && j == 20)
                receiver9(iter_t) = u(i,j);
                iter_t = iter_t + 1;
            end
                        
            % tau2(i,j)  = (tau2(i+1,j+1) + tau2(i+1,j) + tau2(i,j+1) + tau2(i,j)) / 4;
            % M
            % B
            
            % PML for x
            if (i>=nx-nd)
                di=-3*vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)^pow/(2*(nd*dx)^3);
            else
                di=0;
            end
            
            % PML for z
            if (j>=nz-nd)
                dj=-3*vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)^pow/(2*(nd*dz)^3);
            elseif (j<=nd)
                dj=-3*vp(i,j)*log(0.0001)*(nd*dz-j*dz)^pow/(2*(nd*dz)^3);
            else
                dj=0;
            end
            
            % P-SV Wave Propagation, Virieux
            ux(i,j)=(1-dt*di/2)*uxm1(i,j)/(1+dt*di/2)+B(i,j)*dt*(tau1(i+1,j)-tau1(i,j))/((1+dt*di/2)*dx);
            uz(i,j)=(1-dt*dj/2)*uzm1(i,j)/(1+dt*dj/2)+B(i,j)*dt*(tau2(i,j)-tau2(i,j-1))/((1+dt*dj/2)*dz);
            
            vx(i,j)=(1-dt*di/2)*vxm1(i,j)/(1+dt*di/2)+B(i,j)*dt*(tau2(i,j)-tau2(i-1,j))/((1+dt*di/2)*dx);
            vz(i,j)=(1-dt*dj/2)*vzm1(i,j)/(1+dt*dj/2)+B(i,j)*dt*(tau3(i,j+1)-tau3(i,j))/((1+dt*dj/2)*dz);
            
            % u=u(x)+u(z)
            % v=v(x)+v(z)
            
            u(i,j)=ux(i,j)+uz(i,j);
            u(i-1,j)=ux(i-1,j)+uz(i-1,j);
            v(i,j)=vx(i,j)+vz(i,j);
            v(i,j-1)=vx(i,j-1)+vz(i,j-1);
            v(i+1,j)=vx(i+1,j)+vz(i+1,j);
            u(i,j+1)=ux(i,j+1)+uz(i,j+1);
            
            % Stresses
            tau1xp1(i,j)=(1-dt*di/2)*tau1x(i,j)/(1+dt*di/2)+(L(i,j)+2*M(i,j))*dt*(u(i,j)-u(i-1,j))/((1+dt*di/2)*dx);
            tau1zp1(i,j)=(1-dt*dj/2)*tau1z(i,j)/(1+dt*dj/2)+L(i,j)*dt*(v(i,j)-v(i,j-1))/((1+dt*dj/2)*dz);
            
            tau3xp1(i,j)=(1-dt*di/2)*tau3x(i,j)/(1+dt*di/2)+L(i,j)*dt*(u(i,j)-u(i-1,j))/((1+dt*di/2)*dx);
            tau3zp1(i,j)=(1-dt*dj/2)*tau3z(i,j)/(1+dt*dj/2)+(L(i,j)+2*M(i,j))*dt*(v(i,j)-v(i,j-1))/((1+dt*dj/2)*dz);
            
            Mavg = (M(i-1,j) + M(i+1,j) + M(i,j-1) + M(i,j+1))/4;
            tau2xp1(i,j)=(1-dt*di/2)*tau2x(i,j)/(1+dt*di/2)+Mavg*dt*(v(i+1,j)-v(i,j))/((1+dt*di/2)*dx);
            tau2zp1(i,j)=(1-dt*dj/2)*tau2z(i,j)/(1+dt*dj/2)+Mavg*dt*(u(i,j+1)-u(i,j))/((1+dt*dj/2)*dz);
        end
    end
    
    % u(n-1)=u(n)
    uxm1=ux;
    uzm1=uz;
    vxm1=vx;
    vzm1=vz;
    
    % tau(n)=tau(n+1)
    tau1x=tau1xp1;
    tau1z=tau1zp1;
    tau2x=tau2xp1;
    tau2z=tau2zp1;
    tau3x=tau3xp1;
    tau3z=tau3zp1;
    
    % tau=tau(x)+tau(z)
    tau1=tau1xp1+tau1zp1;
    tau2=tau2xp1+tau2zp1;
    tau3=tau3xp1+tau3zp1;
    
    % Reflecting boundary conditions
    %tau1(:,[1 end])=0;
    %tau1([1 end],:)=0;
    %tau3(:,[1 end])=0;
    %tau3([1 end],:)=0;
end

figure();
imagesc(x,z,u);
title(sprintf('Vz; t = %.10f; Nx=%d, Nz=%d, dx=%d, dz=%d, ro=%d', t, nx, nz, dx, dz, ro));
colorbar;
axis equal;
axis([0 Lx 0 Lz]);
xlabel('X'); ylabel('Z');

v=zeros(nx,nz);     % v(n), optional
vx=zeros(nx,nz);    % vx(n)
vz=zeros(nx,nz);    % vz(n)
vxm1=zeros(nx,nz);  % vx(n-1)
vzm1=zeros(nx,nz);  % vz(n-1)
% U
u=zeros(nx,nz);     % u(n), optional
ux=zeros(nx,nz);    % ux(n)
uz=zeros(nx,nz);    % uz(n)
uxm1=zeros(nx,nz);  % ux(n-1)
uzm1=zeros(nx,nz);  % uz(n-1)
% tau(xx)
tau1=zeros(nx,nz);     % tauxx(n), optional
tau1x=zeros(nx,nz);    % tauxx(x)(n)
tau1z=zeros(nx,nz);    % tauxx(z)(n)
tau1xp1=zeros(nx,nz);  % tauxx(x)(n+1)
tau1zp1=zeros(nx,nz);  % tauxx(z)(n+1)
% tau(xz)
tau2=zeros(nx,nz);     % tauxz(n), optional
tau2x=zeros(nx,nz);    % tauxz(x)(n)
tau2z=zeros(nx,nz);    % tauxz(z)(n)
tau2xp1=zeros(nx,nz);  % tauxz(x)(n+1)
tau2zp1=zeros(nx,nz);  % tauxz(z)(n+1)
% tau(zz)
tau3=zeros(nx,nz);     % tauzz(n), optional
tau3x=zeros(nx,nz);    % tauzz(x)(n)
tau3z=zeros(nx,nz);    % tauzz(z)(n)
tau3xp1=zeros(nx,nz);  % tauzz(x)(n+1)
tau3zp1=zeros(nx,nz);  % tauzz(z)(n+1)

t = t - dt;
iter_t = nt;

% Animation
fig = figure();
frame = getframe(fig);
[im,map] = rgb2ind(frame.cdata,4);
imwrite(im,map,'animation.gif','DelayTime',0,'Loopcount',0);

while (t>0)
    t=t-dt;
    % Source
    
    tau3(2,40)=tau3(2,40)+receiver1(iter_t);
    tau3(2,60)=tau3(2,60)+receiver2(iter_t);
    tau3(2,80)=tau3(2,80)+receiver3(iter_t);
    tau3(2,100)=tau3(2,100)+receiver4(iter_t);
    tau3(2,120)=tau3(2,120)+receiver5(iter_t);
    tau3(2,140)=tau3(2,140)+receiver6(iter_t);
    tau3(2,160)=tau3(2,160)+receiver7(iter_t);
    tau3(2,180)=tau3(2,180)+receiver8(iter_t);
    tau3(2,20)=tau3(2,20)+receiver8(iter_t);
    
    iter_t = iter_t - 1;
    
    %tau1(101,101)=tau1(101,101)+source_term*dt;
    %disp(source_term);
    
    for i=2:nx-1
        for j=2:nz-1
            
            % alfaS_out(i,j§)=(alfaS(i,j)+alfaS(i+1,j)+alfaS(i,j+1)+alfaS(i+1,j+1))/4.0d0
            
            % tau2(i,j)  = (tau2(i+1,j+1) + tau2(i+1,j) + tau2(i,j+1) + tau2(i,j)) / 4;
            % M
            % B
            
            % PML for x
            if (i>=nx-nd)
                di=-3*vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)^pow/(2*(nd*dx)^3);
            else
                di=0;
            end
            
            % PML for z
            if (j>=nz-nd)
                dj=-3*vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)^pow/(2*(nd*dz)^3);
            elseif (j<=nd)
                dj=-3*vp(i,j)*log(0.0001)*(nd*dz-j*dz)^pow/(2*(nd*dz)^3);
            else
                dj=0;
            end
            
            % P-SV Wave Propagation, Virieux
            ux(i,j)=(1-dt*di/2)*uxm1(i,j)/(1+dt*di/2)+B(i,j)*dt*(tau1(i+1,j)-tau1(i,j))/((1+dt*di/2)*dx);
            uz(i,j)=(1-dt*dj/2)*uzm1(i,j)/(1+dt*dj/2)+B(i,j)*dt*(tau2(i,j)-tau2(i,j-1))/((1+dt*dj/2)*dz);
            
            vx(i,j)=(1-dt*di/2)*vxm1(i,j)/(1+dt*di/2)+B(i,j)*dt*(tau2(i,j)-tau2(i-1,j))/((1+dt*di/2)*dx);
            vz(i,j)=(1-dt*dj/2)*vzm1(i,j)/(1+dt*dj/2)+B(i,j)*dt*(tau3(i,j+1)-tau3(i,j))/((1+dt*dj/2)*dz);
            
            % u=u(x)+u(z)
            % v=v(x)+v(z)
            
            u(i,j)=ux(i,j)+uz(i,j);
            u(i-1,j)=ux(i-1,j)+uz(i-1,j);
            v(i,j)=vx(i,j)+vz(i,j);
            v(i,j-1)=vx(i,j-1)+vz(i,j-1);
            v(i+1,j)=vx(i+1,j)+vz(i+1,j);
            u(i,j+1)=ux(i,j+1)+uz(i,j+1);
            
            % Stresses
            tau1xp1(i,j)=(1-dt*di/2)*tau1x(i,j)/(1+dt*di/2)+(L(i,j)+2*M(i,j))*dt*(u(i,j)-u(i-1,j))/((1+dt*di/2)*dx);
            tau1zp1(i,j)=(1-dt*dj/2)*tau1z(i,j)/(1+dt*dj/2)+L(i,j)*dt*(v(i,j)-v(i,j-1))/((1+dt*dj/2)*dz);
            
            tau3xp1(i,j)=(1-dt*di/2)*tau3x(i,j)/(1+dt*di/2)+L(i,j)*dt*(u(i,j)-u(i-1,j))/((1+dt*di/2)*dx);
            tau3zp1(i,j)=(1-dt*dj/2)*tau3z(i,j)/(1+dt*dj/2)+(L(i,j)+2*M(i,j))*dt*(v(i,j)-v(i,j-1))/((1+dt*dj/2)*dz);
            
            Mavg = (M(i-1,j) + M(i+1,j) + M(i,j-1) + M(i,j+1))/4;
            tau2xp1(i,j)=(1-dt*di/2)*tau2x(i,j)/(1+dt*di/2)+Mavg*dt*(v(i+1,j)-v(i,j))/((1+dt*di/2)*dx);
            tau2zp1(i,j)=(1-dt*dj/2)*tau2z(i,j)/(1+dt*dj/2)+Mavg*dt*(u(i,j+1)-u(i,j))/((1+dt*dj/2)*dz);
        end
    end
    
    % u(n-1)=u(n)
    uxm1=ux;
    uzm1=uz;
    vxm1=vx;
    vzm1=vz;
    
    % tau(n)=tau(n+1)
    tau1x=tau1xp1;
    tau1z=tau1zp1;
    tau2x=tau2xp1;
    tau2z=tau2zp1;
    tau3x=tau3xp1;
    tau3z=tau3zp1;
    
    % tau=tau(x)+tau(z)
    tau1=tau1xp1+tau1zp1;
    tau2=tau2xp1+tau2zp1;
    tau3=tau3xp1+tau3zp1;
    
    % Reflecting boundary conditions
    %tau1(:,[1 end])=0;
    %tau1([1 end],:)=0;
    %tau3(:,[1 end])=0;
    %tau3([1 end],:)=0;
    
    imagesc(x,z,u);
    title(sprintf('Vz; t = %.10f; Nx=%d, Nz=%d, dx=%d, dz=%d, ro=%d', t, nx, nz, dx, dz, ro));
    colorbar;
    axis equal;
    axis([0 Lx 0 Lz]);
    xlabel('X'); ylabel('Z');
    frame = getframe(fig);
    [im,map] = rgb2ind(frame.cdata,4);
    imwrite(im,map,'animation.gif','DelayTime',0.1,'WriteMode','Append');
end

figure();
imagesc(x,z,u);
title(sprintf('Vz; t = %.10f; Nx=%d, Nz=%d, dx=%d, dz=%d, ro=%d', t, nx, nz, dx, dz, ro));
colorbar;
axis equal;
axis([0 Lx 0 Lz]);
xlabel('X'); ylabel('Z');

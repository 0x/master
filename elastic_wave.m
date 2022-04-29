% clear; close('all')

nx=201;  % Number points of x
nz=201;  % Number points of z

dx=20;  % 20 m
dz=dx;
T=0.3;  % Total time 

vp=zeros(nx,nz);
vs=zeros(nx,nz);

% Read vp clutter model from file
filename = 'Model_clatter_Vp_S_20_I_4';
filenameExtention =[filename '.mod'];
fid = fopen(filenameExtention ,'r');
%vp=fread(fid,[nx,nz],'double');
imagesc(vp);
title('clutter; Nx=Nz=200; hx=hz=20')
colorbar;

ro=2600;
vp(:)=5500;
vs=vp/sqrt(3);  %vs/vp<0.5
%vs=0;

lambda=ro*(vp.^2-2*vs.^2);  % 2.621666666666665e+10
nyu=ro*vs.^2;  % 2.621666666666665e+10

%CFL=0.5;  % CFL=vp*dt/dx 
%dt=CFL*min(dx,dz)/(max(max(vp/2))*sqrt(2.0));

dt = 0.9*dx/(5500*sqrt(2));


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
vm1=zeros(nx,nz);  % vx(n-1)
vzm1=zeros(nx,nz);  % vz(n-1)
% U
u=zeros(nx,nz);     % u(n), optional
ux=zeros(nx,nz);    % ux(n)
uz=zeros(nx,nz);    % uz(n)
um1=zeros(nx,nz);  % ux(n-1)
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
nd=0;  % Thickness of PML
pow = 2;  

% Make constant for now
B(:)=1/ro; 
L(:)=lambda;
M(:)=nyu;

% Scheme
t=0;
while (t<T)
    t=t+dt;
     
    
    
% Source
    source_term=1*(1-2*a*(t-T/2)^2)*exp(-a*(t-T/2)^2); % Rickerra  
    tau3(100,100)=tau3(100,100)+source_term;
    %tau1(100,100)=tau1(100,100)+source_term*dt;
    %disp(source_term);
    
    
    for i=2:nx-1
        for j=2:nz-1           

            % P-SV Wave Propagation, Virieux
            u(i,j)=u(i,j)+B(i,j)*dt*(tau1(i+1,j)-tau1(i,j))/dx ...
                +B(i,j)*dt*(tau2(i,j)-tau2(i,j-1))/dx;
          
            v(i,j)=v(i,j)+B(i,j)*dt*(tau2(i,j)-tau2(i-1,j))/dx ...
                +B(i,j)*dt*(tau3(i,j+1)-tau3(i,j))/dx;
                    % Stresses
            tau1xp1(i,j)=tau1(i,j)+(L(i,j)+2*M(i,j))*dt*(u(i,j)-u(i-1,j))/dx ...
                +L(i,j)*dt*(v(i,j)-v(i,j-1))/dx;
           
            tau3xp1(i,j)=tau3(i,j)+(L(i,j)+2*M(i,j))*dt*(v(i,j)-v(i,j-1))/dx ... 
                +L(i,j)*dt*(u(i,j)-u(i-1,j))/dx;
            
            Mavg = (M(i-1,j) + M(i+1,j) + M(i,j-1) + M(i,j+1))/4;
            tau2xp1(i,j)=tau2(i,j)+Mavg*dt*(v(i+1,j)-v(i,j))/dx ...
                +Mavg*dt*(u(i,j+1)-u(i,j))/dx;
        end
    end
    
    % tau(n)=tau(n+1)
    tau1=tau1xp1;
    tau2=tau2xp1;
    tau3=tau3xp1;
     u([1 nx],:)   =  0;
    u(:,[1 nz])   =  0;
    v([1 nx],:)   =  0;
    v(:,[1 nz])   =  0;  
    % Reflecting boundary conditions
    tau1(:,[1 end])=0;
    tau1([1 end],:)=0;
    tau3(:,[1 end])=0;
    tau3([1 end],:)=0;
end

figure();
imagesc(x,z,v); 
title(sprintf('Vz; t = %.10f; Nx=%d, Nz=%d, dx=%d, dz=%d, ro=%d', t, nx, nz, dx, dz, ro));
colorbar;
axis equal; 
axis([0 Lx 0 Lz]); 
xlabel('X'); ylabel('Z');

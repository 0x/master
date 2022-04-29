program elastic_wave

implicit none

! Source in the global coordinate
integer, parameter :: isrc = 101, jsrc = 101

! Wave equation parameters
real(8), parameter :: dx = 20, dz = 20
real(8), parameter :: vp = 5500, ro = 2600, t = 0.3, vs = 0
real(8), parameter :: t0 = t / 2
real(8), parameter :: CFL = 0.5
real(8), parameter :: dt = CFL * min(dx, dz) / ((vp / 2) * sqrt(2.0))
real(8) :: tt = 0
integer :: i, j

! Wave source
real(8), parameter :: pi = 3.141592653589793238462643d0, f0 = 15
real(8), parameter :: a = pi * pi * f0 * f0
real(8) :: source_term

! Lame coefficients
real(8), parameter :: lambda = ro * (vp * vp - 2 * vs * vs), nyu = ro * vs * vs

! Domain
integer, parameter :: nx = 201, nz = 201
real(8), parameter :: Lx = nx * dx, Ly = nz * dz

! PML
real(8) :: di = 0, dj = 0
integer, parameter :: nd = 20  ! Thickness of PML
integer, parameter :: pow = 2

real(8), dimension(:,:), allocatable :: v, vx, vz, vxm1, vzm1, u, ux, uz, uxm1, uzm1
real(8), dimension(:,:), allocatable :: tau1, tau1x, tau1z, tau1xp1, tau1zp1
real(8), dimension(:,:), allocatable :: tau2, tau2x, tau2z, tau2xp1, tau2zp1
real(8), dimension(:,:), allocatable :: tau3, tau3x, tau3z, tau3xp1, tau3zp1
real(8), dimension(:,:), allocatable :: B, L, M

! Uz
allocate(v(nx, nz));
allocate(vx(nx, nz));
allocate(vz(nx, nz));
allocate(vxm1(nx, nz));
allocate(vzm1(nx, nz));
! Ux
allocate(u(nx, nz));
allocate(ux(nx, nz));
allocate(uz(nx, nz));
allocate(uxm1(nx, nz));
allocate(uzm1(nx, nz));
! tauxx
allocate(tau1(nx, nz));
allocate(tau1x(nx, nz));
allocate(tau1z(nx, nz));
allocate(tau1xp1(nx, nz));
allocate(tau1zp1(nx, nz));
! tauxz
allocate(tau2(nx, nz));
allocate(tau2x(nx, nz));
allocate(tau2z(nx, nz));
allocate(tau2xp1(nx, nz));
allocate(tau2zp1(nx, nz));
! tauzz
allocate(tau3(nx, nz));
allocate(tau3x(nx, nz));
allocate(tau3z(nx, nz));
allocate(tau3xp1(nx, nz));
allocate(tau3zp1(nx, nz));

allocate(B(nx, nz));
allocate(L(nx, nz));
allocate(M(nx, nz));

! Init
v = 0
vx = 0
vz = 0
vxm1 = 0
vzm1 = 0
u = 0
ux = 0
uz = 0
uxm1 = 0
uzm1 = 0
tau1 = 0
tau1x = 0
tau1z = 0
tau1xp1 = 0
tau1zp1 = 0
tau2 = 0
tau2x = 0
tau2z = 0
tau2xp1 = 0
tau2zp1 = 0
tau3 = 0
tau3x = 0
tau3z = 0
tau3xp1 = 0
tau3zp1 = 0

! Const for now
B = 1 / ro
L = lambda
M = nyu

!Calculation of the wave equation
do while (tt < t)
	tt = tt + dt

	uxm1 = ux
	uzm1 = uz
	vxm1 = vx
	vzm1 = vz

	tau1x = tau1xp1
	tau1z = tau1zp1
	tau2x = tau2xp1
	tau2z = tau2zp1
	tau3x = tau3xp1
	tau3z = tau3zp1

	tau1 = tau1x + tau1z
	tau2 = tau2x + tau2z
	tau3 = tau3x + tau3z

	source_term = 1 * (1.0 - 2.0 * a * (tt - t0)**2)*exp(-a * (tt - t0)**2)
	tau3(101, 101) = tau3(101, 101) + source_term * dt

	do i = 2,nx-1
		do j = 2,nz-1
			! PML for x
			if (i<=nd) then
				di=-3*vp*log(0.0001)*(nd*dx-i*dx)**pow/(2*(nd*dx)**3)
			else if (i>=nx-nd) then
				di=-3*vp*log(0.0001)*((+i-nx+nd)*dx)**pow/(2*(nd*dx)**3)
			else
				di=0
			end if
			! PML for z
			if (j>=nz-nd) then
				dj=-3*vp*log(0.0001)*((+j-nz+nd)*dz)**pow/(2*(nd*dz)**3)
			else if (j<=nd) then
				dj=-3*vp*log(0.0001)*(nd*dz-j*dz)**pow/(2*(nd*dz)**3)
			else
				dj=0
			end if
	 
			ux(i,j)=(1-dt*di/2)*uxm1(i,j)/(1+dt*di/2)+B(i,j)*dt*(tau1(i+1,j)-tau1(i,j))/((1+dt*di/2)*dx)
	 		uz(i,j)=(1-dt*dj/2)*uzm1(i,j)/(1+dt*dj/2)+B(i,j)*dt*(tau2(i,j)-tau2(i,j-1))/((1+dt*dj/2)*dz)
	 
	 		vx(i,j)=(1-dt*di/2)*vxm1(i,j)/(1+dt*di/2)+B(i,j)*dt*(tau2(i,j)-tau2(i-1,j))/((1+dt*di/2)*dx)
	 		vz(i,j)=(1-dt*dj/2)*vzm1(i,j)/(1+dt*dj/2)+B(i,j)*dt*(tau3(i,j+1)-tau3(i,j))/((1+dt*dj/2)*dz)
	 
	 		u(i,j)=ux(i,j)+uz(i,j)
	 		u(i-1,j)=ux(i-1,j)+uz(i-1,j)
	 		v(i,j)=vx(i,j)+vz(i,j)
	 		v(i,j-1)=vx(i,j-1)+vz(i,j-1)
	 		v(i+1,j)=vx(i+1,j)+vz(i+1,j)
	 		u(i,j+1)=ux(i,j+1)+uz(i,j+1)

	 		tau1xp1(i,j)=(1-dt*di/2)*tau1x(i,j)/(1+dt*di/2)+(L(i,j)+2*M(i,j))*dt*(u(i,j)-u(i-1,j))/((1+dt*di/2)*dx)
	 		tau1zp1(i,j)=(1-dt*dj/2)*tau1z(i,j)/(1+dt*dj/2)+L(i,j)*dt*(v(i,j)-v(i,j-1))/((1+dt*dj/2)*dz)
	  
	 		tau3xp1(i,j)=(1-dt*di/2)*tau3x(i,j)/(1+dt*di/2)+L(i,j)*dt*(u(i,j)-u(i-1,j))/((1+dt*di/2)*dx)
	 		tau3zp1(i,j)=(1-dt*dj/2)*tau3z(i,j)/(1+dt*dj/2)+(L(i,j)+2*M(i,j))*dt*(v(i,j)-v(i,j-1))/((1+dt*dj/2)*dz)
   
	 		tau2xp1(i,j)=(1-dt*di/2)*tau2x(i,j)/(1+dt*di/2)+M(i,j)*dt*(v(i+1,j)-v(i,j))/((1+dt*di/2)*dx)
	 		tau2zp1(i,j)=(1-dt*dj/2)*tau2z(i,j)/(1+dt*dj/2)+M(i,j)*dt*(u(i,j+1)-u(i,j))/((1+dt*dj/2)*dz)
		enddo
	enddo
enddo

open(unit = 1 , file = "elastic_wave.txt")
write(1,*) v
close(1)
	
end program elastic_wave

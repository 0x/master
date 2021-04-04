!mpirun --mca btl_vader_backing_directory /tmp --oversubscribe  -np 4 elastic_mpi

program elastic_wave

use mpi

implicit none

! MPI stuff
integer :: error, nprocs, rank
integer :: istatus(MPI_STATUS_SIZE)
integer, dimension(MPI_STATUS_SIZE) :: status

! DD
integer, parameter :: iprocs = 2, jprocs = 2

! Source in the global coordinate
integer, parameter :: isrc = 50, jsrc = 50
integer :: local_isrc, local_jsrc

! Wave equation parameters
real(8), parameter :: dx = 20.0, dz = 20.0
real(8), parameter :: ro = 2000.0, t = 1.0
real(8), parameter :: t0 = t / 2.0
real(8), parameter :: CFL = 0.5
real(8) :: dt
real(8) :: tt = 0.0
integer :: i, j

! Wave source
real(8), parameter :: pi = 3.141592653589793238462643d0, f0 = 15.0
real(8), parameter :: a = pi * pi * f0 * f0
real(8) :: source_term

! Domain
integer, parameter :: nxlocal = 100, nzlocal = 100
integer, parameter :: nx=nxlocal+2, nz=nzlocal+2
real(8), parameter :: Lx=nxlocal*dx, Lz=nzlocal*dz

integer, parameter :: sizeOfReal8 = 8
real(8), dimension(:,:), allocatable :: vp, vs, lambda, nyu

! PML
real(8) :: di = 0.0, dj = 0.0
integer, parameter :: nd = 20  ! Thickness of PML
integer, parameter :: pow = 2

real(8), dimension(:,:), allocatable :: v, vx, vz, vxm1, vzm1, u, ux, uz, uxm1, uzm1
real(8), dimension(:,:), allocatable :: tau1, tau1x, tau1z, tau1xp1, tau1zp1
real(8), dimension(:,:), allocatable :: tau2, tau2x, tau2z, tau2xp1, tau2zp1
real(8), dimension(:,:), allocatable :: tau3, tau3x, tau3z, tau3xp1, tau3zp1
real(8), dimension(:,:), allocatable :: B, L, M
real(8) :: Mavg
integer :: ista, iend, jsta, jend

! Itable interface
integer :: irank, ranki, rankj
integer :: inext, jnext, iprev, jprev

real(8), dimension(:), allocatable :: tau2_s, tau1_s, vx_s, vz_s, ux_s, uz_s
real(8), dimension(:), allocatable :: tau2_r, tau1_r, vx_r, vz_r, ux_r, uz_r
real(8), dimension(:), allocatable :: send1, send2, send3, send4
real(8), dimension(:), allocatable :: recv1, recv2, recv3, recv4

! MPI I/O
integer(kind=MPI_OFFSET_KIND), parameter :: zero_off = 0
integer :: globalsize(2), subsize(2), starts(2), locallenght, written_arr, read_arr, fhandle_w, fhandle_r

character*60 :: input_Vp_filename="Model_clatter_Vp_S_20_I_4.mod"
character*60 :: outfile
integer :: index
integer :: iter_snap=0, num_snap=15
character*60, external :: index_of_file

! MPI init
call mpi_init(error)
call mpi_comm_size(MPI_COMM_WORLD, nprocs, error)
call mpi_comm_rank(MPI_COMM_WORLD, rank, error)

if (nprocs /= iprocs*jprocs) then
    write(*, *) 'ERROR, mpi_comm_size: ', nprocs
    stop
end if

call calculate_itable(iprocs, jprocs, rank, inext, jnext, iprev, jprev, ranki, rankj) 
    
! TODO: decrease allocations 
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

allocate(vp(nx, nz));
allocate(vs(nx, nz));
allocate(lambda(nx, nz));
allocate(nyu(nx, nz));

allocate(send1(3*nx)); allocate(send2(3*nx)); allocate(recv1(3*nx)); allocate(recv2(3*nx));
allocate(send3(3*nz)); allocate(send4(3*nz)); allocate(recv3(3*nz)); allocate(recv4(3*nz));

! Initialization
v = 0.0
vx = 0.0
vz = 0.0
vxm1 = 0.0
vzm1 = 0.0
u = 0.0
ux = 0.0
uz = 0.0
uxm1 = 0.0
uzm1 = 0.0
tau1 = 0.0
tau1x = 0.0
tau1z = 0.0
tau1xp1 = 0.0
tau1zp1 = 0.0
tau2 = 0.0
tau2x = 0.0
tau2z = 0.0
tau2xp1 = 0.0
tau2zp1 = 0.0
tau3 = 0.0
tau3x = 0.0
tau3z = 0.0
tau3xp1 = 0.0
tau3zp1 = 0.0
vp = 0.0
vs = 0.0
send1 = 0.0
send2 = 0.0
send3 = 0.0
send4 = 0.0
recv1 = 0.0
recv2 = 0.0
recv3 = 0.0
recv4 = 0.0

! MPI I/O
globalsize(1)=nxlocal*iprocs
globalsize(2)=nzlocal*jprocs
subsize(1)=nxlocal
subsize(2)=nzlocal
starts(1)= ranki*nxlocal
starts(2)= rankj*nzlocal
locallenght=subsize(1)*subsize(2)
call mpi_type_create_subarray(2, globalsize, subsize, starts, MPI_ORDER_FORTRAN, MPI_REAL8, written_arr, error)
call mpi_type_commit(written_arr, error)
call mpi_type_create_subarray(2, globalsize, subsize, starts, MPI_ORDER_FORTRAN, MPI_REAL8, read_arr, error)
call mpi_type_commit(read_arr, error)

call mpi_file_open(MPI_COMM_WORLD, input_Vp_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fhandle_r, error)
call mpi_file_set_view(fhandle_r, zero_off, MPI_REAL8, read_arr,  "native", MPI_INFO_NULL, error)
call mpi_file_read_all(fhandle_r, vp(2:nx-1,2:nz-1), locallenght, MPI_REAL8, istatus, error)
call mpi_barrier(MPI_COMM_WORLD, error)
call mpi_file_close(fhandle_r, error)

! Calculate Vs from Vp 
vs = 0!vp/sqrt(3.0)

! TODO: make dt for all ranks
! Calculate time step
!dt=CFL*min(dx,dz)/((maxval(vp/2.0))*sqrt(2.0))
dt = 0.0046
write(*, *) dt

! Lame coeficients
lambda = ro * (vp * vp - 2.0 * vs * vs)
nyu = ro * vs * vs

B = 1.0 / ro
L = lambda
M = nyu

! Calculation of the wave equation
do while (tt < t)
    tt = tt + dt

    send1(1:nx) = tau3(:,2)
    send1(nx+1:2*nx) = ux(:,2)
    send1(2*nx+1:3*nx) = uz(:,2)
    
    send2(1:nx) = tau2(:,nz-1)
    send2(nx+1:2*nx) = vx(:,nz-1)
    send2(2*nx+1:3*nx) = vz(:,nz-1)
    
    send3(1:nz) = tau1(2,:)
    send3(nz+1:2*nz) = vx(2,:)
    send3(2*nz+1:3*nz) = vz(2,:)

    send4(1:nz) = tau2(nx-1,:)
    send4(nz+1:2*nz) = ux(nx-1,:)
    send4(2*nz+1:3*nz) = uz(nx-1,:)

    call mpi_send(send2, 3*nx, MPI_REAL8, jnext, 1, MPI_COMM_WORLD, error)
    call mpi_send(send1, 3*nx, MPI_REAL8, jprev, 1, MPI_COMM_WORLD, error)
        
    call mpi_send(send4, 3*nz, MPI_REAL8, inext, 1, MPI_COMM_WORLD, error)
    call mpi_send(send3, 3*nz, MPI_REAL8, iprev, 1, MPI_COMM_WORLD, error)
    
    call mpi_recv(recv2, 3*nx, MPI_REAL8, jprev, 1, MPI_COMM_WORLD, status, error)
    call mpi_recv(recv1, 3*nx, MPI_REAL8, jnext, 1, MPI_COMM_WORLD, status, error)

    call mpi_recv(recv4, 3*nz, MPI_REAL8, iprev, 1, MPI_COMM_WORLD, status, error)
    call mpi_recv(recv3, 3*nz, MPI_REAL8, inext, 1, MPI_COMM_WORLD, status, error)

    tau3(:,nz) = recv1(1:nx)
    ux(:,nz) = recv1(nx+1:2*nx)
    uz(:,nz) = recv1(2*nx+1:3*nx)

    tau2(:,1) = recv2(1:nx)
    vx(:,1) = recv2(nx+1:2*nx)
    vz(:,1) = recv2(2*nx+1:3*nx)

    tau1(nx,:) = recv3(1:nz)
    vx(nx,:) = recv3(nz+1:2*nz)
    vz(nx,:) = recv3(2*nz+1:3*nz)

    tau2(1,:) = recv4(1:nz)
    ux(1,:) = recv4(nz+1:2*nz)
    uz(1,:) = recv4(2*nz+1:3*nz)

    ! TODO: make subroutine
    ! source 
    if (isrc .gt. ranki*nxlocal .and. isrc .lt. (ranki+1)*nxlocal .and. &
        jsrc .gt. rankj*nzlocal .and. jsrc .lt. (rankj+1)*nzlocal) then 
        source_term = 1.0 * (1.0 - 2.0 * a * (tt - t0)**2)*exp(-a * (tt - t0)**2)
        local_isrc = isrc-ranki*nxlocal
        local_jsrc = jsrc-rankj*nzlocal
        !write(*,*) local_isrc, local_jsrc
        tau3(local_isrc, local_jsrc) = tau3(local_isrc, local_jsrc) + source_term * dt
    end if 

    ista=2
    jsta=2
    iend=nx-1
    jend=nz-1
    if (ranki == 0) then
        ista=ista+1
    endif
    if (ranki == iprocs - 1) then
        iend=iend-1
    endif
    if (rankj == 0) then
        jsta=jsta+1
    endif
    if (rankj == jprocs - 1) then
        jend=jend-1
    endif

    do i = ista,iend
        do j = jsta,jend
            ! PML for x
            if (ranki*nxlocal+i >= iprocs*nxlocal-nd ) then
                di=-3*vp(i,j)*log(0.0001)*((i-nx+nd)*dx)**pow/(2*(nd*dx)**3)
            else
                di=0
            end if
            ! PML for z
            if (rankj*nzlocal+j>=jprocs*nzlocal-nd) then
                dj=-3*vp(i,j)*log(0.0001)*((j-nz+nd)*dz)**pow/(2*(nd*dz)**3)
            else if (rankj*nzlocal+j<=nd) then
                dj=-3*vp(i,j)*log(0.0001)*(nd*dz-j*dz)**pow/(2*(nd*dz)**3)
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
            
            ! Calculate M(i+1/2, j+1/2) as average 
            Mavg=(M(i-1,j)+M(i+1,j)+M(i,j-1)+M(i,j+1))/4.0; 

            tau2xp1(i,j)=(1-dt*di/2)*tau2x(i,j)/(1+dt*di/2)+Mavg*dt*(v(i+1,j)-v(i,j))/((1+dt*di/2)*dx)
            tau2zp1(i,j)=(1-dt*dj/2)*tau2z(i,j)/(1+dt*dj/2)+Mavg*dt*(u(i,j+1)-u(i,j))/((1+dt*dj/2)*dz)
        end do
    end do

    ! V(n-1) = V(n)
    uxm1 = ux
    uzm1 = uz
    vxm1 = vx
    vzm1 = vz

    ! tau(n) = tau(n+1)
    tau1x = tau1xp1
    tau1z = tau1zp1
    tau2x = tau2xp1
    tau2z = tau2zp1
    tau3x = tau3xp1
    tau3z = tau3zp1

    tau1 = tau1xp1 + tau1zp1
    tau2 = tau2xp1 + tau2zp1
    tau3 = tau3xp1 + tau3zp1
end do

call mpi_file_open(MPI_COMM_WORLD, "test_mpi.mod", MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fhandle_w, error)
call mpi_file_set_view(fhandle_w, zero_off, MPI_REAL8, written_arr,  "native", MPI_INFO_NULL, error)
call mpi_file_write_all(fhandle_w, u(2:nx-1,2:nz-1), locallenght, MPI_REAL8, istatus, error)
call mpi_barrier(MPI_COMM_WORLD, error)
call mpi_file_close(fhandle_w, error)

call mpi_finalize(error)

contains

subroutine calculate_itable(iprocs, jprocs, rank, inext, jnext, iprev, jprev, myranki, myrankj)     
    
    integer, intent(in) :: iprocs, jprocs, rank
    integer, intent(out) :: inext, jnext, iprev, jprev, myranki, myrankj

    integer :: itable(-1:iprocs, -1:jprocs)
    integer :: irank, i, j

    do j = -1, jprocs
        do i = -1, iprocs
            itable(i,j) = MPI_PROC_NULL
        end do
    end do
    irank = 0
    do i = 0, iprocs-1
        do j = 0, jprocs-1
            itable(i,j) = irank
            if (rank == irank) then
                myranki = i;
                myrankj = j
            end if
            irank = irank + 1
        end do
    end do
    jnext = itable(myranki, myrankj + 1)
    jprev = itable(myranki, myrankj - 1)
    inext = itable(myranki+1, myrankj)
    iprev = itable(myranki-1, myrankj)

end subroutine calculate_itable
    
end program elastic_wave

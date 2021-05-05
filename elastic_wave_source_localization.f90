program elastic_wave

use mpi

implicit none

! MPI
integer :: error, nprocs, rank
integer :: istatus(MPI_STATUS_SIZE)
integer, dimension(MPI_STATUS_SIZE) :: status

integer, parameter :: iprocs = 4, jprocs = 4

! Source in the global coordinate
integer, parameter :: isrc = 51, jsrc = 101

! Wave equation parameters
real(8), parameter :: dx = 20.0, dz = 20.0
real(8), parameter :: r = 2600.0, time_global = 1.3
real(8), parameter :: t0 = time_global / 2
real(8), parameter :: CFL = 0.5
real(8) :: dt, a
real(8) :: time = 0.0
integer :: i, j, nt, n

! Wave source
real(8), parameter :: pi = 3.141592653589793238462643d0, f0 = 15.0
real(8), parameter :: ac = pi * pi * f0 * f0
real(8) :: source_term

! Domain
integer, parameter :: nx = 50, nz = 50
integer, parameter :: nx_mpi = nx+2, nz_mpi = nz+2
integer, parameter :: sizeOfReal8 = 8
real(8), dimension(:,:), allocatable :: Vp, Vs, lambda, miu, rho
real(8), parameter :: Lx = nx * dx, Ly = nz * dz

! PML
real(8) :: di = 0.0, dj = 0.0
integer, parameter :: nd = 20  ! Thickness of PML
integer, parameter :: pow = 2

real(8), dimension(:,:), allocatable :: v, vx, vz, u, ux, uz
real(8), dimension(:,:), allocatable :: Txx, Txx_x, Txx_z, Tzz, Tzz_x, Tzz_z, Txz, Txz_x, Txz_z
real(8), dimension(:,:), allocatable :: BU, LAM, MU, GAMMA, b


! Itable interface
integer :: irank, ranki, rankj
integer :: inext, jnext, iprev, jprev
integer :: ista, iend, jsta, jend
integer :: local_isrc, local_jsrc

real(8), dimension(:), allocatable :: Txz_send, Txz_recv, Txx_send, Txx_recv
real(8), dimension(:), allocatable :: u_send, u_recv, v_send, v_recv

real(8), dimension(:), allocatable :: receiver1, receiver2, receiver3, receiver4, receiver5
real(8), dimension(:), allocatable :: receiver6, receiver7, receiver8, receiver9
real(8) :: epsilon_xx, epsilon_xz, epsilon_zz, nu, E
real(8), dimension(:,:), allocatable :: energy, max_energy

! MPI I/O
integer(kind=MPI_OFFSET_KIND), parameter :: zero_off = 0
integer :: globalsize(2), subsize(2), starts(2), locallenght, written_arr, read_arr, fhandle_w, fhandle_r

character*60 :: input_Vp_filename="Model_clatter_Vp_S_20_I_4.mod"
character*60 :: outfile
integer :: index
integer :: iter_snap=0, num_snap=15
character*60, external :: index_of_file 

call mpi_init(error)
call mpi_comm_size(MPI_COMM_WORLD, nprocs, error)
call mpi_comm_rank(MPI_COMM_WORLD, rank, error)

call calculate_itable(iprocs, jprocs, rank, inext, jnext, iprev, jprev, ranki, rankj) 

! Read Vp from bin file
globalsize(1)=nx*iprocs
globalsize(2)=nz*jprocs
subsize(1)=nx
subsize(2)=nz
starts(1)= ranki*nx
starts(2)= rankj*nz
locallenght=subsize(1)*subsize(2)
call mpi_type_create_subarray(2, globalsize, subsize, starts, MPI_ORDER_FORTRAN, MPI_REAL8, written_arr, error)
call mpi_type_commit(written_arr, error)
call mpi_type_create_subarray(2, globalsize, subsize, starts, MPI_ORDER_FORTRAN, MPI_REAL8, read_arr, error)
call mpi_type_commit(read_arr, error)

allocate(Vp(nx, nz))
call mpi_file_open(MPI_COMM_WORLD, input_Vp_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fhandle_r, error)
call mpi_file_set_view(fhandle_r, zero_off, MPI_REAL8, read_arr,  "native", MPI_INFO_NULL, error)
call mpi_file_read_all(fhandle_r, vp, locallenght, MPI_REAL8, istatus, error)
call mpi_barrier(MPI_COMM_WORLD, error)
call mpi_file_close(fhandle_r, error)

! Uz
allocate(Vs(nx, nz))
Vs = Vp/sqrt(3.0)

allocate(rho(nx, nz))
allocate(miu(nx, nz))
allocate(lambda(nx, nz))
rho = r
miu = (Vs*Vs)*r
lambda = (Vp*Vp)*r-2.0*miu

! TODO: make equal dt to all ranks
dt = 0.9*dx/(5500.0*sqrt(2.0));

allocate(u(0:nx+1, 0:nz+1)) 
allocate(ux(nx, nz))
allocate(uz(nx, nz))
allocate(v(0:nx+1, 0:nz+1))
allocate(vx(nx, nz))
allocate(vz(nx, nz))

allocate(Txx(0:nx+1, 0:nz+1))
allocate(Tzz(0:nx+1, 0:nz+1))
allocate(Txz(0:nx+1, 0:nz+1))
allocate(Txx_x(nx, nz)); allocate(Txx_z(nx, nz))
allocate(Tzz_x(nx, nz)); allocate(Tzz_z(nx, nz))
allocate(Txz_x(nx, nz)); allocate(Txz_z(nx, nz))

u = 0.0
ux = 0.0
uz = 0.0
v = 0.0
vx = 0.0
vz = 0.0

Txx = 0.0
Tzz = 0.0
Txz = 0.0
Txx_x = 0.0
Tzz_x = 0.0
Txz_x = 0.0
Txx_z = 0.0
Tzz_z = 0.0
Txz_z = 0.0

! Bouyancy and other parameter
allocate(b(nx, nz))
b = 1.0/rho;
a = dt/dx;

allocate(BU(nx, nz))
allocate(LAM(nx, nz))
allocate(MU(nx, nz))
allocate(GAMMA(nx, nz))
BU = b*a;
LAM = lambda*a;
MU = miu*a;
GAMMA = LAM + 2.0*MU;

! Stations
nt = int(time_global/dt)+1
allocate(receiver1(nt))
allocate(receiver2(nt))
allocate(receiver3(nt))
allocate(receiver4(nt))
allocate(receiver5(nt))
allocate(receiver6(nt))
allocate(receiver7(nt))
allocate(receiver8(nt))
allocate(receiver9(nt))
allocate(energy(nx, nz))
allocate(max_energy(nx, nz))


allocate(Txz_send(nz_mpi)); allocate(Txz_recv(nz_mpi))
allocate(Txx_send(nz_mpi)); allocate(Txx_recv(nz_mpi))
allocate(u_send(nz_mpi)); allocate(u_recv(nz_mpi))
allocate(v_send(nz_mpi)); allocate(v_recv(nz_mpi))

! Calculation of the wave equation
n = 1
do while (time < time_global)

    !Txz(:,1) = Txz(:,nz_mpi-1)
    !Tzz(:,nz_mpi) = Tzz(:,2)
    !Txx(nx+1,:)=Txx(1,:)
    !Txz(0,:)=Txz(nx,:)

    Txz_send = Txz(nx,:)
    Txx_send = Txx(1,:)

    call mpi_send(Txz(:,nz), nx_mpi, MPI_REAL8, jnext, 1000, MPI_COMM_WORLD, error)
    call mpi_send(Tzz(:,1), nx_mpi, MPI_REAL8, jprev, 2000, MPI_COMM_WORLD, error)
    
    call mpi_send(Txz_send, nz_mpi, MPI_REAL8, inext, 5000, MPI_COMM_WORLD, error)
    call mpi_send(Txx_send, nz_mpi, MPI_REAL8, iprev, 6000, MPI_COMM_WORLD, error)

    call mpi_recv(Txz_recv, nz_mpi, MPI_REAL8, iprev, 5000, MPI_COMM_WORLD, status, error)
    call mpi_recv(Txx_recv, nz_mpi, MPI_REAL8, inext, 6000, MPI_COMM_WORLD, status, error)

    call mpi_recv(Txz(:,0), nx_mpi, MPI_REAL8, jprev, 1000, MPI_COMM_WORLD, status, error)
    call mpi_recv(Tzz(:,nz+1), nx_mpi, MPI_REAL8, jnext, 2000, MPI_COMM_WORLD, status, error)
    
    if (ranki /= 0) then
        Txz(0,:) = Txz_recv
    end if
    if (ranki /= iprocs-1) then
        Txx(nx+1,:) = Txx_recv
    endif

	do i=1,nx
        do j=1,nz
            ! PML for x
            if (ranki*nx+i >= iprocs*nx-nd) then
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)**pow/(2*(nd*dx)**3)
            else
                di=0
            endif
            ! PML for z
            if (rankj*nz+j>=jprocs*nz-nd) then
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)**pow/(2*(nd*dz)**3)
            else if (rankj*nz+j<=nd) then
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)**pow/(2*(nd*dz)**3)
            else
                dj=0
            endif
            
            ux(i,j) = (1-dt*di/2)*ux(i,j)/(1+dt*di/2) + BU(i,j) * (Txx(i+1,j)-Txx(i,j))/(1+dt*di/2)
            uz(i,j) = (1-dt*dj/2)*uz(i,j)/(1+dt*dj/2) + BU(i,j) * (Txz(i,j)-Txz(i,j-1))/(1+dt*dj/2)
            u(i,j) = ux(i,j) + uz(i,j)
            
            vx(i,j) = (1-dt*di/2)*vx(i,j)/(1+dt*di/2) + BU(i,j) * (Txz(i,j)-Txz(i-1,j))/(1+dt*di/2)
            vz(i,j) = (1-dt*dj/2)*vz(i,j)/(1+dt*dj/2) + BU(i,j) * (Tzz(i,j+1)- Tzz(i,j))/(1+dt*dj/2)
            v(i,j) = vx(i,j) + vz(i,j)

            ! TODO: make uniform distribution
            if (i == 2 .and. ranki == 0) then 
                select case (j+rankj*nz)
                case (40)
                    receiver1(n) = u(i,j)
                case (60)
                    receiver2(n) = u(i,j)
                case (80)
                    receiver3(n) = u(i,j)
                case (100)
                    receiver4(n) = u(i,j)
                case (120)
                    receiver5(n) = u(i,j)
                case (140)
                    receiver6(n) = u(i,j)
                case (160)
                    receiver7(n) = u(i,j)
                case (180)
                    receiver8(n) = u(i,j)
                case (20)
                    receiver9(n) = u(i,j)
                end select
            endif
        
        enddo
    enddo

    !v(:,1) = v(:,nz_mpi-1)
    !u(:,nz_mpi) = u(:,2)
    !u(0,:)=u(nx,:)
    !v(nx+1,:)=v(1,:)

    u_send = u(nx,:)
    v_send = v(1,:)

    call mpi_send(v(:,nz), nx_mpi, MPI_REAL8, jnext, 3000, MPI_COMM_WORLD, error)
    call mpi_send(u(:,1), nx_mpi, MPI_REAL8, jprev, 4000, MPI_COMM_WORLD, error)

    call mpi_send(u_send, nz_mpi, MPI_REAL8, inext, 7000, MPI_COMM_WORLD, error)
    call mpi_send(v_send, nz_mpi, MPI_REAL8, iprev, 8000, MPI_COMM_WORLD, error)
    
    call mpi_recv(v(:,0), nx_mpi, MPI_REAL8, jprev, 3000, MPI_COMM_WORLD, status, error)
    call mpi_recv(u(:,nz+1), nx_mpi, MPI_REAL8, jnext, 4000, MPI_COMM_WORLD, status, error)

    call mpi_recv(u_recv, nz_mpi, MPI_REAL8, iprev, 7000, MPI_COMM_WORLD, status, error)
    call mpi_recv(v_recv, nz_mpi, MPI_REAL8, inext, 8000, MPI_COMM_WORLD, status, error)
 
    if (ranki /= 0) then
        u(0,:) = u_recv
    end if
    if (ranki /= iprocs-1) then
        v(nx+1,:) = v_recv
    endif

    if (ranki*nx < isrc .and. isrc <= (ranki+1)*nx .and. rankj*nz < jsrc .and. jsrc <= (rankj+1)*nz) then
        source_term=1*(1-2*ac*(time-t0)**2)*exp(-ac*(time-t0)**2); ! Rickerra
        
        Tzz_z(mod(isrc, ranki*nx),mod(jsrc, rankj*nz)) =  Tzz_z(mod(isrc, ranki*nx),mod(jsrc, rankj*nz))+source_term*dt;
        !Txx_x(mod(isrc, ranki*nx),mod(jsrc, rankj*nz)) =  Txx_x(mod(isrc, ranki*nx),mod(jsrc, rankj*nz))+source_term*dt;
    endif

    do i=1,nx
        do j=1,nz
            ! PML for x
            if (ranki*nx+i >= iprocs*nx-nd) then
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)**pow/(2*(nd*dx)**3)
            else
                di=0
            endif
            ! PML for z
            if (rankj*nz+j>=jprocs*nz-nd) then
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)**pow/(2*(nd*dz)**3)
            else if (rankj*nz+j<=nd) then
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)**pow/(2*(nd*dz)**3)
            else
                dj=0
            endif
            
            Txx_x(i,j) = (1-dt*di/2)*Txx_x(i,j)/(1+dt*di/2) + GAMMA(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2)
            Txx_z(i,j) = (1-dt*dj/2)*Txx_z(i,j)/(1+dt*dj/2) + LAM(i,j)*(v(i,j)-v(i,j-1))/(1+dt*dj/2)
            Txx(i,j) = Txx_x(i,j) + Txx_z(i,j)
            
            Tzz_x(i,j) = (1-dt*di/2)*Tzz_x(i,j)/(1+dt*di/2) + LAM(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2)
            Tzz_z(i,j) = (1-dt*dj/2)*Tzz_z(i,j)/(1+dt*dj/2) + GAMMA(i,j)*(v(i,j)- v(i,j-1))/(1+dt*dj/2)
            Tzz(i,j) = Tzz_x(i,j) + Tzz_z(i,j)
            
            Txz_x(i,j) = (1-dt*di/2)*Txz_x(i,j)/(1+dt*di/2) + MU(i,j) * (v(i+1,j)- v(i,j))/(1+dt*di/2)
            Txz_z(i,j) = (1-dt*dj/2)*Txz_z(i,j)/(1+dt*dj/2) + MU(i,j)*(u(i,j+1)-u(i,j))/(1+dt*dj/2)
            Txz(i,j) = Txz_x(i,j) + Txz_z(i,j)
        enddo
    enddo

    time = time+dt
    n = n + 1
enddo

call mpi_file_open(MPI_COMM_WORLD, "wave.mod", MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fhandle_w, error)
call mpi_file_set_view(fhandle_w, zero_off, MPI_REAL8, written_arr,  "native", MPI_INFO_NULL, error)
call mpi_file_write_all(fhandle_w, v(1:nx,1:nz), locallenght, MPI_REAL8, istatus, error)
call mpi_barrier(MPI_COMM_WORLD, error)
call mpi_file_close(fhandle_w, error)

! Calculation of the wave equation with receivers as source
u = 0.0
ux = 0.0
uz = 0.0
v = 0.0
vx = 0.0
vz = 0.0

Txx = 0.0
Tzz = 0.0
Txz = 0.0
Txx_x = 0.0
Tzz_x = 0.0
Txz_x = 0.0
Txx_z = 0.0
Tzz_z = 0.0
Txz_z = 0.0

n = nt
time = 0.0
do while (time < time_global)

    !Txz(:,1) = Txz(:,nz_mpi-1)
    !Tzz(:,nz_mpi) = Tzz(:,2)
    !Txx(nx+1,:)=Txx(1,:)
    !Txz(0,:)=Txz(nx,:)

    Txz_send = Txz(nx,:)
    Txx_send = Txx(1,:)

    call mpi_send(Txz(:,nz), nx_mpi, MPI_REAL8, jnext, 1000, MPI_COMM_WORLD, error)
    call mpi_send(Tzz(:,1), nx_mpi, MPI_REAL8, jprev, 2000, MPI_COMM_WORLD, error)
    
    call mpi_send(Txz_send, nz_mpi, MPI_REAL8, inext, 5000, MPI_COMM_WORLD, error)
    call mpi_send(Txx_send, nz_mpi, MPI_REAL8, iprev, 6000, MPI_COMM_WORLD, error)

    call mpi_recv(Txz_recv, nz_mpi, MPI_REAL8, iprev, 5000, MPI_COMM_WORLD, status, error)
    call mpi_recv(Txx_recv, nz_mpi, MPI_REAL8, inext, 6000, MPI_COMM_WORLD, status, error)

    call mpi_recv(Txz(:,0), nx_mpi, MPI_REAL8, jprev, 1000, MPI_COMM_WORLD, status, error)
    call mpi_recv(Tzz(:,nz+1), nx_mpi, MPI_REAL8, jnext, 2000, MPI_COMM_WORLD, status, error)
    
    if (ranki /= 0) then
        Txz(0,:) = Txz_recv
    end if
    if (ranki /= iprocs-1) then
        Txx(nx+1,:) = Txx_recv
    endif

    do i=1,nx
        do j=1,nz
            ! PML for x
            if (ranki*nx+i >= iprocs*nx-nd) then
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)**pow/(2*(nd*dx)**3)
            else
                di=0
            endif
            ! PML for z
            if (rankj*nz+j>=jprocs*nz-nd) then
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)**pow/(2*(nd*dz)**3)
            else if (rankj*nz+j<=nd) then
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)**pow/(2*(nd*dz)**3)
            else
                dj=0
            endif
            
            ux(i,j) = (1-dt*di/2)*ux(i,j)/(1+dt*di/2) + BU(i,j) * (Txx(i+1,j)-Txx(i,j))/(1+dt*di/2)
            uz(i,j) = (1-dt*dj/2)*uz(i,j)/(1+dt*dj/2) + BU(i,j) * (Txz(i,j)-Txz(i,j-1))/(1+dt*dj/2)
            u(i,j) = ux(i,j) + uz(i,j)
            
            vx(i,j) = (1-dt*di/2)*vx(i,j)/(1+dt*di/2) + BU(i,j) * (Txz(i,j)-Txz(i-1,j))/(1+dt*di/2)
            vz(i,j) = (1-dt*dj/2)*vz(i,j)/(1+dt*dj/2) + BU(i,j) * (Tzz(i,j+1)- Tzz(i,j))/(1+dt*dj/2)
            v(i,j) = vx(i,j) + vz(i,j)
        enddo
    enddo

    !v(:,1) = v(:,nz_mpi-1)
    !u(:,nz_mpi) = u(:,2)
    !u(0,:)=u(nx,:)
    !v(nx+1,:)=v(1,:)

    u_send = u(nx,:)
    v_send = v(1,:)

    call mpi_send(v(:,nz), nx_mpi, MPI_REAL8, jnext, 3000, MPI_COMM_WORLD, error)
    call mpi_send(u(:,1), nx_mpi, MPI_REAL8, jprev, 4000, MPI_COMM_WORLD, error)

    call mpi_send(u_send, nz_mpi, MPI_REAL8, inext, 7000, MPI_COMM_WORLD, error)
    call mpi_send(v_send, nz_mpi, MPI_REAL8, iprev, 8000, MPI_COMM_WORLD, error)
    
    call mpi_recv(v(:,0), nx_mpi, MPI_REAL8, jprev, 3000, MPI_COMM_WORLD, status, error)
    call mpi_recv(u(:,nz+1), nx_mpi, MPI_REAL8, jnext, 4000, MPI_COMM_WORLD, status, error)

    call mpi_recv(u_recv, nz_mpi, MPI_REAL8, iprev, 7000, MPI_COMM_WORLD, status, error)
    call mpi_recv(v_recv, nz_mpi, MPI_REAL8, inext, 8000, MPI_COMM_WORLD, status, error)
 
    if (ranki /= 0) then
        u(0,:) = u_recv
    end if
    if (ranki /= iprocs-1) then
        v(nx+1,:) = v_recv
    endif
    
    ! TODO: rewrite with subroutines
    if (ranki == 0) then
        if (rankj*nz < 40 .and.  40 <= (rankj+1)*nz) &
            Tzz_z(2,40)=Tzz_z(2,40)+receiver1(n)
        if (rankj*nz < 60 .and. 60 <= (rankj+1)*nz) &
            Tzz_z(2,mod(60, rankj*nz))=Tzz_z(2,mod(60, rankj*nz))+receiver3(n)
        if (rankj*nz < 80 .and. 80 <= (rankj+1)*nz) &
            Tzz_z(2,mod(80, rankj*nz))=Tzz_z(2,mod(80, rankj*nz))+receiver3(n)
        if (rankj*nz < 100 .and. 100 <= (rankj+1)*nz) &
            Tzz_z(2,mod(100, rankj*nz))=Tzz_z(2,mod(100, rankj*nz))+receiver4(n)
        if (rankj*nz < 120 .and. 120 <= (rankj+1)*nz) &
            Tzz_z(2,mod(120, rankj*nz))=Tzz_z(2,mod(120, rankj*nz))+receiver5(n)
        if (rankj*nz < 140 .and. 140 <= (rankj+1)*nz) &
            Tzz_z(2,mod(140, rankj*nz))=Tzz_z(2,mod(140, rankj*nz))+receiver6(n)
        if (rankj*nz < 160 .and. 160 <= (rankj+1)*nz) &
            Tzz_z(2,mod(160, rankj*nz))=Tzz_z(2,mod(160, rankj*nz))+receiver7(n)
        if (rankj*nz < 180 .and. 180 <= (rankj+1)*nz) &
            Tzz_z(2,mod(180, rankj*nz))=Tzz_z(2,mod(180, rankj*nz))+receiver8(n)
        if (rankj*nz < 20 .and. 20 <= (rankj+1)*nz) &
            Tzz_z(2,20)=Tzz_z(2,20)+receiver9(n)
    endif

    do i=1,nx
        do j=1,nz
            ! PML for x
            if (ranki*nx+i >= iprocs*nx-nd) then
                di=-3*Vp(i,j)*log(0.0001)*((+i-nx+nd)*dx)**pow/(2*(nd*dx)**3)
            else
                di=0
            endif
            ! PML for z
            if (rankj*nz+j>=jprocs*nz-nd) then
                dj=-3*Vp(i,j)*log(0.0001)*((+j-nz+nd)*dz)**pow/(2*(nd*dz)**3)
            else if (rankj*nz+j<=nd) then
                dj=-3*Vp(i,j)*log(0.0001)*(nd*dz-j*dz)**pow/(2*(nd*dz)**3)
            else
                dj=0
            endif
            
            Txx_x(i,j) = (1-dt*di/2)*Txx_x(i,j)/(1+dt*di/2) + GAMMA(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2)
            Txx_z(i,j) = (1-dt*dj/2)*Txx_z(i,j)/(1+dt*dj/2) + LAM(i,j)*(v(i,j)-v(i,j-1))/(1+dt*dj/2)
            Txx(i,j) = Txx_x(i,j) + Txx_z(i,j)
            
            Tzz_x(i,j) = (1-dt*di/2)*Tzz_x(i,j)/(1+dt*di/2) + LAM(i,j)*(u(i,j)-u(i-1,j))/(1+dt*di/2)
            Tzz_z(i,j) = (1-dt*dj/2)*Tzz_z(i,j)/(1+dt*dj/2) + GAMMA(i,j)*(v(i,j)- v(i,j-1))/(1+dt*dj/2)
            Tzz(i,j) = Tzz_x(i,j) + Tzz_z(i,j)
            
            Txz_x(i,j) = (1-dt*di/2)*Txz_x(i,j)/(1+dt*di/2) + MU(i,j) * (v(i+1,j)- v(i,j))/(1+dt*di/2)
            Txz_z(i,j) = (1-dt*dj/2)*Txz_z(i,j)/(1+dt*dj/2) + MU(i,j)*(u(i,j+1)-u(i,j))/(1+dt*dj/2)
            Txz(i,j) = Txz_x(i,j) + Txz_z(i,j)
        enddo
    enddo

    ! Calculate energy
    do i=1,nx
        do j=1,nz
            nu = 0.5 * lambda(i,j)/(lambda(i,j) + miu(i,j));
            E = 2.0*miu(i,j)*(1.0+ nu);
            ! Compute total field from split components
            epsilon_xx = ( Txx(i,j) - nu*Tzz(i,j) )/E;
            epsilon_zz = ( Tzz(i,j) - nu*Txx(i,j) )/E;
            epsilon_xz = 2.0*(1.0+ nu)*Txz(i,j) / E;
            
            energy(i,j) = 0.50 * ( 0.5*(r+r)* u(i,j)**2 + 0.5*(r+r)*v(i,j)**2) + &
                0.50 * ( epsilon_xx * Txx(i,j) + epsilon_zz * Tzz(i,j) + epsilon_xz * Txz(i,j));
            
            max_energy(i, j) = (max_energy(i, j)+ energy(i, j));
        enddo
    enddo

    time = time+dt
    n = n - 1
enddo

call mpi_file_open(MPI_COMM_WORLD, "max_energy.mod", MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fhandle_w, error)
call mpi_file_set_view(fhandle_w, zero_off, MPI_REAL8, written_arr,  "native", MPI_INFO_NULL, error)
call mpi_file_write_all(fhandle_w, max_energy(1:nx,1:nz), locallenght, MPI_REAL8, istatus, error)
call mpi_barrier(MPI_COMM_WORLD, error)
call mpi_file_close(fhandle_w, error)

call mpi_file_open(MPI_COMM_WORLD, "v.mod", MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fhandle_w, error)
call mpi_file_set_view(fhandle_w, zero_off, MPI_REAL8, written_arr,  "native", MPI_INFO_NULL, error)
call mpi_file_write_all(fhandle_w, v(1:nx,1:nz), locallenght, MPI_REAL8, istatus, error)
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

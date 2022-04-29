   program hello
   use mpi
   integer rank, size, ierror, tag, status(MPI_STATUS_SIZE)
   
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
   print*, 'node', rank, ': Hello world'
   call MPI_FINALIZE(ierror)
   end

   if (ranki == 0) then
        if (rankj*nz < 40 .and.  40 <= (rankj+1)*nz) &
            Tzz_z(2,40)=Tzz_z(2,40)+receiver1(n)
        if (rankj*nz < 60 .and. 60 <= (rankj+1)*nz) &
            write(*,*) mod(60, rankj*nz)
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
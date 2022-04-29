program elastic_wave

implicit none

integer, parameter :: nx=200, nz=200
real(8), dimension(nx, nz) :: vp


open(unit=13, file='Model_clatter_Vp_S_20_I_4.mod', form='unformatted', access='direct', recl=nx*nz*8)       
read(13, rec=1) vp   
close (13)

write(*,*) vp

end program elastic_wave

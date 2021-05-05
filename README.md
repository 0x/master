# master

## Example
- Setup mpifort elastic_wave_clatter_mpi.f90

Set processes number in -L 12 ```integer, parameter :: iprocs = 4, jprocs = 4``` 
- Build

```mpifort elastic_wave_clatter_mpi.f90 -o elastic_wave_clatter_mpi```
- Run 

```mpirun --mca btl_vader_backing_directory /tmp --oversubscribe  -np 16 ./elastic_wave_clatter_mpi```

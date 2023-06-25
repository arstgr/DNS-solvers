# Quadrant analysis
This repository contains source codes for performing quadrant analysis in turbulent channel flows. The computations are performed using the statistical analysis suggested by [Wallace et al. 1972](https://doi.org/10.1017/S0022112072000515) which has become a standard tool in the study of [wall turbulence since then](https://doi.org/10.1146/annurev-fluid-122414-034550).  

The program is designed to perform computations in parallel. It is essential that the number of grid points in the streamwise direction is divisible by the total number of MPI ranks.

Parallel I/O operations are employed by the program. All input and output files are expected to reside on a shared file system.

## Getting Started

### Dependencies

* None.

### To Build
  Simply type
```
  make 
```

### To Run
The code is designed to run in a parallel computing environment. Use the following command to execute the program:

* The code runs in parallel, using mpi:
```
mpirun -np N qudranta
```
where N denotes the number of MPI ranks.

### Input files
The program expects the following input files:

  * "vel.time" (e.g. vel.0001) : velocity file containing the instantaneous velocity field at time T, with array indices running from (0..Nx-1, 0..Ny-1, 0..Nz-1, 0..2), where the indices go by first z index, then y index, then x, then i (velocity vector component). The files are expected to be row major, as is the case for the C language.

  * "den.time" (e.g. den.0001) : density file containing the instantaneous density (pressure) field at time T, with array indices running from (0..Nx-1, 0..Ny-1, 0..Nz-1). 

  * The input parameters for the computations also need to be modified at the top of `src/definitions.h` file.

#### Output files
The program generates the following output file:

quadrant.000.time corresponding to quadrant analysis results (i.e. 4 components) at every z plane, first for the velocity then for the pressure, at the given time T.

#### Postprocessing

To perform statistical analysis on the resulting data sets, try typing
```
make postprocess
```
then run the resulting binary 
```
./postprocess
```

The postprocessing component requires the pressure gradient and flow rate at every time step of the simulation, to calculate the wall shear stress and the resulting friction velocity and also the mean bulk flow velocity. The output *.dat files can be opened in most general purpose visualization software, including tecplot(TM) and gnuplot. 

## Author

Amirreza Rastegari [@arstgr](https://github.com/arstgr)

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the GNU Affero general public License. Refer to the [LICENSE](LICENSE) file for more information.


## To Cite

If you use this code in your research, please cite the following work:

Rastegari, S.A., 2017. Computational studies of turbulent skin-friction drag reduction with super-hydrophobic surfaces and riblets ([Doctoral dissertation](https://deepblue.lib.umich.edu/handle/2027.42/136986)), the University of Michigan.
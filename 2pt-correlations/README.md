# 2-Point Velocity Auto-Correlations
This repository contains source codes for calculating the 2-Point velocity auto-correlations in turbulent channel flows. 

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
mpirun -np N 2ptcorl
```
where N denotes the number of MPI ranks.

### Input files
The program expects the following input file:

"vel.time" (e.g. vel.0001) : velocity file containing the instantaneous velocity field at time T, with array indices running from (0..Nx-1, 0..Ny-1, 0..Nz-1), where the indices go by first z index, then y index, then x 

#### Output files
The program generates the following output files:

xcorl.time and ycorl.time corresponding to 2-Point velocity auto-correlations in the streamwise (x) and spanwise (y) directions, respectively, at the given time T.

#### Postprocessing

To perform statistical analysis on the resulting data sets, first modify the input parameters at the top of `src/postproc-corl.c`, then build the postprocessing component by typing
```
make postprocess
```
then run the resulting binary 
```
./postprocess
```

The postprocessing component requires the pressure gradient at every time step of the simulation, to calculate the wall shear stress and the resulting friction velocity. Additionally, the beginning and final time during which the statistical averages are to be calculated, along with the channel sizes need to be specified at the top of the source file. The output *.dat files can be opened in most general purpose visualization software, including tecplot(TM) and gnuplot. 

## Results

Due to the smaller channel sizes in my DNS studies, the correlations specially the one in the spanwise direction, don't go perfectly to zero. However they match those published by Jimenez et al. and Moser et al. which were performed at much larger channel sizes. The spanwise auto correlations near the wall accurately predict the mean size of streamwise vortices as well, as seen in the figures bellow. The Figures demonstrate 2-Point velocity auto-correlations obtained, using this code, at a distance of 10 wall units from the walls, in turbulent channel flow at a bulk Reynolds number of Reb=Ubxh/nu=7200. 
![2-Point velocity auto-correlations in streamwise and spanwise directions](figs/2pt-corrls.png)


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
# Turbulence statistics and budget of turbulence kinetic energy

This repository contains the source codes for computing the turbulence statistics and budget of turbulence kinetic energy in turbulent channel flow with micro-textured walls. The computations are done in parallel, while the postprocessing is done in serial.

## Requirements

- MPI 
- C compiler 

## Installation

1. Clone the repository and CD to the repo
2. Modify the grid sizes and problem specifications in the `definitions.h` file located in the respective source directories:
   - `computations/definitions.h`: Modify the grid sizes and problem parameters for the computations.
   - `postprocessing/definitions.h`: Modify the problem parameters for postprocessing.

3. Compile the codes by running the following command:
```
    make
```

## Usage

To perform the computations:
```
mpirun -np N ./compute-budget
```

Replace `N` with the desired number of MPI processes for parallel computation. This command will execute the compute program, which will perform the computations using 4th order or 2nd order compact finite difference stencils. The results will be stored in output files.

To perform the postprocessing:
```
./postprocess_budget
```
This command will execute the postprocess program, which will process the results obtained from the computations. The postprocessing is performed in serial.

### Input files
The program expects the following input file:

  * "vel.time" (e.g. vel.0001) : velocity file containing the instantaneous velocity field at time T, with array indices running from (0..Nx-1, 0..Ny-1, 0..Nz-1, 0..2), where the indices go by first z index, then y index, then x, then i (velocity vector component). The files are expected to be row major, as is the case for the C language.

  * "den.time" (e.g. den.0001) : density file containing the instantaneous density (pressure) field at time T, with array indices running from (0..Nx-1, 0..Ny-1, 0..Nz-1). 

  * "p-grad.txt" time history of the evolution of the pressure gradient applied to the flow (either constant or time dependent) 

  * "flowrate.txt" time history of the flow rate in the channel (either constant or time dependent)

#### Output files
The program outputs are formatted files with ".dat" extension that can be read into general purpose visualization software. 

## Results

The computations all use 4th order or second order compact finite difference stencils everywhere in the domain, and hence the overall accuracy of the results is second order. Results obtained from these codes have been published in the journal of fluid mechanics, [here](https://doi.org/10.1017/jfm.2015.266) and [here](https://doi.org/10.1017/jfm.2017.865). 


The results are in full agreement with the published DNS databases of Jimenez et al. and Moser et al

## Author

Amirreza Rastegari [@arstgr](https://github.com/arstgr)

## To Cite

If you use this code in your research, please cite the following work:

Rastegari, S.A., 2017. Computational studies of turbulent skin-friction drag reduction with super-hydrophobic surfaces and riblets ([Doctoral dissertation](https://deepblue.lib.umich.edu/handle/2027.42/136986)), the University of Michigan.

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the GNU Affero general public License. Refer to the [LICENSE](LICENSE) file for more information.

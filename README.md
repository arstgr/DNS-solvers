# DNS-solvers
Parallel solvers for Direct Numerical Simulation (DNS) of turbulent channel flow and postprocessing of the results

This repository contains a set of parallel solvers for the Direct Numerical Simulation (DNS) of turbulent channel flows. The solvers are all second order in space and time, and provide excellent scalability on thousands of cores. The solvers allow any periodic pattern of micro-texture and indentations to be placed on the walls. The repository is structured as follows

  * [DNS-Solver](https://github.com/arstgr/DNS-solvers/tree/main/DNS-solver)
  * [TKE-budget](https://github.com/arstgr/DNS-solvers/tree/main/TKE-budget)
  * [1D Energy spectra](https://github.com/arstgr/DNS-solvers/tree/main/e-spectra)
  * [2-Point velocity auto correlations](https://github.com/arstgr/DNS-solvers/tree/main/2pt-correlations)
  * [Quadrant analysis](https://github.com/arstgr/DNS-solvers/tree/main/quadrant-analysis)

![Contour plot of Coherent Structures](figs/fig1.jpg) 

## To Cite

If you use this code in your research, please cite the following work:

  * Rastegari, S.A., 2017. Computational studies of turbulent skin-friction drag reduction with super-hydrophobic surfaces and riblets ([Doctoral dissertation](https://deepblue.lib.umich.edu/handle/2027.42/136986)), the University of Michigan.

  * Rastegari, A. and Akhavan, R., 2015. On the mechanism of turbulent drag reduction with super-hydrophobic surfaces, JFM, 773 (R4). 

  * Rastegari, A. and Akhavan, R., 2018. The Common Mechanism of Turbulent Skin-Friction Drag Reduction with Super-Hydrophobic Micro-Grooves and Riblets, JFM, 838 (68-104). 

# advanced-computational-engineering-fortran
📌 **Overview:**

This repository contains a collection of numerical simulation projects implemented in Fortran, focusing on heat transfer, diffusion problems, and high-performance computing (HPC) using MPI and OpenMP.
The work demonstrates the application of finite difference methods, iterative solvers, and parallel computing techniques to solve engineering problems efficiently.

📂 **Projects:**
**1. 2D Bi-material Heat Diffusion (MPI + OpenMP):**

Simulation of heat transfer in a two-dimensional composite plate, solved using finite difference methods and iterative schemes.

* Bi-material domain with varying diffusivity
* Time-dependent heat diffusion
* Implemented using:
  * Serial code
  * OpenMP parallelisation
  * Hybrid MPI + OpenMP

📄 **Report:** [View Full Report](docs/2d-bimaterial-heat-diffusion/2d-bimaterial-heat-diffusion-report.pdf)
**2. 2D Advection-Diffusion Solver (MPI):**

Parallel solution of a 2D steady advection-diffusion equation using:

* Finite difference discretisation
* Jacobi / Red-Black iterative solver
* Domain decomposition with MPI
* Parallel efficiency analysis

📄 **Report:** [View Full Report](docs/2d-heat-diffusion_mpi/2d-heat-diffusion_mpi-report.pdf)

🛠️ **Tools & Technologies:**
* Fortran
* MPI (OpenMPI)
* OpenMP
* Finite Difference Methods (FDM)
* Iterative Solvers (Jacobi, Red-Black)
  
📁 **Repository Structure**
src/        # Source code (Fortran)
docs/       # Reports

🎯 **Key Highlights:**
* Development of parallel solvers using MPI and OpenMP
* Implementation of finite difference methods for PDEs
* Analysis of parallel efficiency and scalability
* Application to heat transfer and diffusion problems
___________________________________________________________________________________________________________________________________________________________________

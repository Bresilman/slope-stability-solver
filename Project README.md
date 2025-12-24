Slope Stability Solver V1.7

A high-performance Fortran-based engineering tool for slope stability analysis using Limit Equilibrium Methods (LEM).

Features

Methods: Fellenius (Ordinary) and Bishop Simplified.

Reinforcement: Supports Soil Nails (Type 2) and Geogrids (Type 1) with customizable pullout and bond strengths.

Surcharges: Multiple area loads can be applied to the slope surface.

Search: Global minimum search using a radius-center grid approach.

Hydrology: Supports water table definition and pore pressure calculation.

Prerequisites

gfortran (GNU Fortran compiler)

make (Build automation tool)

How to Build

To compile the program, navigate to the project root and run:

make


How to Run

Ensure your geometry and soil properties are defined in input.txt, then execute:

./slope_solver


Project Structure

src/: Contains .f90 source modules and the main program.

build/: Target directory for compiled .o and .mod files.

input.txt: Main configuration file.

results.txt: Summary output of the critical Factor of Safety.
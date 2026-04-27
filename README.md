# 1D Frame Analysis Solver (FEM)

A Python-based implementation of a solver for 1D frame structures with circular cross-sections using the Finite Element Method (FEM).

## Project Overview
The solver follows a standard FEM workflow to:
- Parse structural input data (nodes, elements, materials, boundary conditions, and loads)
- Compute element length and orientation
- Form local stiffness matrices and transform them to the global coordinate system
- Assemble the global system of equations and solve for nodal displacements
- Compute internal forces, strain, and stress, along with visualization of the deformed structure

## My Role
As part of a collaborative team, I focused on understanding and documenting the computational workflow. I studied how the global stiffness matrix is assembled and how the numerical solution is obtained, and helped organize the code and explanation of the overall pipeline.

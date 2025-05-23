# FEMT - Finite Element Method Toolkit

## Introduction

FEMT (Finite Element Method Toolkit) is a computational mechanics finite element analysis system that includes pre-processing, computation, and post-processing components. The original project only supported Intel Fortran compiler. This version has been modified to support:

1. 9-node quadrilateral elements (RECTANGLE9)
2. GNU Fortran (gfortran) compiler

## Project Structure

- `Source/`: Source code directory (computation program)
- `pre_process.m`: Pre-processing program (MATLAB)
- `post_process.m`: Post-processing program (MATLAB)
- `compile.bat`: gfortran compilation script

## Modifications

See `modification.md` for details about the changes made to support 9-node elements and gfortran compilation.

## Usage Instructions

1. Pre-processing: Run `pre_process.m` in MATLAB to generate mesh and boundary conditions
2. Computation: Use `compile.bat` to compile the program, then run the generated executable
3. Post-processing: Run `post_process.m` in MATLAB to visualize results

## Compilation

This project uses GNU Fortran (gfortran) compiler and has been tested with:
```
gcc version 6.3.0 (MinGW.org GCC-6.3.0-1)
```

To compile, use the provided batch file:

```bash
compile.bat
```

The compiled executable will be located in the `gfortran_build` directory.

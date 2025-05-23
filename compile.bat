@echo off
echo Compiling FEMT with gfortran...

REM Create new folders for gfortran build
echo Creating directories for gfortran build...
if not exist "gfortran_build" mkdir gfortran_build
if not exist "gfortran_build\obj" mkdir gfortran_build\obj
if not exist "gfortran_build\obj\Modules" mkdir gfortran_build\obj\Modules
if not exist "gfortran_build\obj\Main" mkdir gfortran_build\obj\Main
if not exist "gfortran_build\obj\Main\Initiation" mkdir gfortran_build\obj\Main\Initiation
if not exist "gfortran_build\obj\Main\Read" mkdir gfortran_build\obj\Main\Read
if not exist "gfortran_build\obj\Boundary" mkdir gfortran_build\obj\Boundary
if not exist "gfortran_build\obj\Force" mkdir gfortran_build\obj\Force
if not exist "gfortran_build\obj\Mechanics" mkdir gfortran_build\obj\Mechanics
if not exist "gfortran_build\obj\Mechanics\Stiffness" mkdir gfortran_build\obj\Mechanics\Stiffness
if not exist "gfortran_build\obj\Mechanics\Strain" mkdir gfortran_build\obj\Mechanics\Strain
if not exist "gfortran_build\obj\Mechanics\Stress" mkdir gfortran_build\obj\Mechanics\Stress
if not exist "gfortran_build\obj\Solver" mkdir gfortran_build\obj\Solver

REM Set compilation options
set FFLAGS=-O2 -J"gfortran_build\obj\Modules" -I"gfortran_build\obj\Modules" -fbacktrace -ffree-line-length-none -fno-range-check -std=legacy

REM First compile module files
echo Compiling modules...
gfortran %FFLAGS% -c Source\Modules\BASIC_DATA.F90 -o gfortran_build\obj\Modules\BASIC_DATA.o
gfortran %FFLAGS% -c Source\Modules\SOLUTION_DATA.F90 -o gfortran_build\obj\Modules\SOLUTION_DATA.o

REM Compile files in Main directory
echo Compiling Main files...
gfortran %FFLAGS% -c Source\Main\CLEAR_INSTRUCTION.F90 -o gfortran_build\obj\Main\CLEAR_INSTRUCTION.o
gfortran %FFLAGS% -c Source\Main\OUTPUT.F90 -o gfortran_build\obj\Main\OUTPUT.o
gfortran %FFLAGS% -c Source\Main\SOLVER_MANAGER.F90 -o gfortran_build\obj\Main\SOLVER_MANAGER.o
gfortran %FFLAGS% -c Source\Main\SUMMARY.F90 -o gfortran_build\obj\Main\SUMMARY.o

REM Compile files in Main\Initiation directory
echo Compiling Initiation files...
gfortran %FFLAGS% -c Source\Main\Initiation\CHECK_DATA.F90 -o gfortran_build\obj\Main\Initiation\CHECK_DATA.o
gfortran %FFLAGS% -c Source\Main\Initiation\INIT_ELE_LIB.F90 -o gfortran_build\obj\Main\Initiation\INIT_ELE_LIB.o
gfortran %FFLAGS% -c Source\Main\Initiation\INIT_INTEGRATION.F90 -o gfortran_build\obj\Main\Initiation\INIT_INTEGRATION.o
gfortran %FFLAGS% -c Source\Main\Initiation\INIT_MECHANICAL.F90 -o gfortran_build\obj\Main\Initiation\INIT_MECHANICAL.o
gfortran %FFLAGS% -c Source\Main\Initiation\INIT_SHAPE.F90 -o gfortran_build\obj\Main\Initiation\INIT_SHAPE.o
gfortran %FFLAGS% -c Source\Main\Initiation\INIT_SOLUTION.F90 -o gfortran_build\obj\Main\Initiation\INIT_SOLUTION.o
gfortran %FFLAGS% -c Source\Main\Initiation\INITIATE.F90 -o gfortran_build\obj\Main\Initiation\INITIATE.o
gfortran %FFLAGS% -c Source\Main\Initiation\LOAD_SHAPE.F90 -o gfortran_build\obj\Main\Initiation\LOAD_SHAPE.o
gfortran %FFLAGS% -c Source\Main\Initiation\SHAPE_2D.F90 -o gfortran_build\obj\Main\Initiation\SHAPE_2D.o

REM Compile files in Main\Read directory
echo Compiling Read files...
gfortran %FFLAGS% -c Source\Main\Read\READ_ELEMENTS.F90 -o gfortran_build\obj\Main\Read\READ_ELEMENTS.o
gfortran %FFLAGS% -c Source\Main\Read\READ_GEOMETRIES.F90 -o gfortran_build\obj\Main\Read\READ_GEOMETRIES.o
gfortran %FFLAGS% -c Source\Main\Read\READ_INSTRUCTION.F90 -o gfortran_build\obj\Main\Read\READ_INSTRUCTION.o
gfortran %FFLAGS% -c Source\Main\Read\READ_MATERIALS.F90 -o gfortran_build\obj\Main\Read\READ_MATERIALS.o
gfortran %FFLAGS% -c Source\Main\Read\READ_XYZ.F90 -o gfortran_build\obj\Main\Read\READ_XYZ.o

REM Compile files in Boundary directory
echo Compiling Boundary files...
gfortran %FFLAGS% -c Source\Boundary\BOUNDARY_CONDITION.F90 -o gfortran_build\obj\Boundary\BOUNDARY_CONDITION.o
gfortran %FFLAGS% -c Source\Boundary\BOUNDARY_DISPLACEMENT.F90 -o gfortran_build\obj\Boundary\BOUNDARY_DISPLACEMENT.o
gfortran %FFLAGS% -c Source\Boundary\BOUNDARY_FORCE.F90 -o gfortran_build\obj\Boundary\BOUNDARY_FORCE.o
gfortran %FFLAGS% -c Source\Boundary\READ_DISPLACEMENT.F90 -o gfortran_build\obj\Boundary\READ_DISPLACEMENT.o
gfortran %FFLAGS% -c Source\Boundary\READ_LINE_FORCE.F90 -o gfortran_build\obj\Boundary\READ_LINE_FORCE.o
gfortran %FFLAGS% -c Source\Boundary\READ_NODAL_FORCE.F90 -o gfortran_build\obj\Boundary\READ_NODAL_FORCE.o

REM Compile files in Force directory
echo Compiling Force files...
gfortran %FFLAGS% -c Source\Force\GET_FORCE.F90 -o gfortran_build\obj\Force\GET_FORCE.o
gfortran %FFLAGS% -c Source\Force\LINE_FORCE_PLANE.F90 -o gfortran_build\obj\Force\LINE_FORCE_PLANE.o

REM Compile files in Mechanics\Stiffness directory
echo Compiling Stiffness files...
gfortran %FFLAGS% -c Source\Mechanics\Stiffness\MATERIAL_2D.F90 -o gfortran_build\obj\Mechanics\Stiffness\MATERIAL_2D.o
gfortran %FFLAGS% -c Source\Mechanics\Stiffness\STIFFNESS.F90 -o gfortran_build\obj\Mechanics\Stiffness\STIFFNESS.o
gfortran %FFLAGS% -c Source\Mechanics\Stiffness\STIFFNESS_PLANE.F90 -o gfortran_build\obj\Mechanics\Stiffness\STIFFNESS_PLANE.o

REM Compile files in Mechanics\Strain directory
echo Compiling Strain files...
gfortran %FFLAGS% -c Source\Mechanics\Strain\GET_STRAIN.F90 -o gfortran_build\obj\Mechanics\Strain\GET_STRAIN.o
gfortran %FFLAGS% -c Source\Mechanics\Strain\STRAIN_GAUSS_PLANE.F90 -o gfortran_build\obj\Mechanics\Strain\STRAIN_GAUSS_PLANE.o
gfortran %FFLAGS% -c Source\Mechanics\Strain\STRAIN_NODAL_RECTANGLE.F90 -o gfortran_build\obj\Mechanics\Strain\STRAIN_NODAL_RECTANGLE.o
gfortran %FFLAGS% -c Source\Mechanics\Strain\STRAIN_NODAL_TRIANGLE.F90 -o gfortran_build\obj\Mechanics\Strain\STRAIN_NODAL_TRIANGLE.o

REM Compile files in Mechanics\Stress directory
echo Compiling Stress files...
gfortran %FFLAGS% -c Source\Mechanics\Stress\GET_STRESS.F90 -o gfortran_build\obj\Mechanics\Stress\GET_STRESS.o
gfortran %FFLAGS% -c Source\Mechanics\Stress\STRESS_GAUSS_PLANE.F90 -o gfortran_build\obj\Mechanics\Stress\STRESS_GAUSS_PLANE.o
gfortran %FFLAGS% -c Source\Mechanics\Stress\STRESS_NODAL_RECTANGLE.F90 -o gfortran_build\obj\Mechanics\Stress\STRESS_NODAL_RECTANGLE.o
gfortran %FFLAGS% -c Source\Mechanics\Stress\STRESS_NODAL_TRIANGLE.F90 -o gfortran_build\obj\Mechanics\Stress\STRESS_NODAL_TRIANGLE.o

REM Compile files in Solver directory
echo Compiling Solver files...
gfortran %FFLAGS% -c Source\Solver\ASSEMBLE_STIFFNESS.F90 -o gfortran_build\obj\Solver\ASSEMBLE_STIFFNESS.o
gfortran %FFLAGS% -c Source\Solver\CHOLESCKEY.F90 -o gfortran_build\obj\Solver\CHOLESCKEY.o
gfortran %FFLAGS% -c Source\Solver\INV.f90 -o gfortran_build\obj\Solver\INV.o
gfortran %FFLAGS% -c Source\Solver\SOLVE.F90 -o gfortran_build\obj\Solver\SOLVE.o
gfortran %FFLAGS% -c Source\Solver\SOLVER_MECH_STATIC_DIR.F90 -o gfortran_build\obj\Solver\SOLVER_MECH_STATIC_DIR.o

REM Compile main program
echo Compiling main program...
gfortran %FFLAGS% -c Source\Main\FEMT.F90 -o gfortran_build\obj\Main\FEMT.o

REM Link all object files to generate executable
echo Linking...
gfortran -o gfortran_build\FEMT.exe ^
gfortran_build\obj\Modules\BASIC_DATA.o ^
gfortran_build\obj\Modules\SOLUTION_DATA.o ^
gfortran_build\obj\Main\CLEAR_INSTRUCTION.o ^
gfortran_build\obj\Main\OUTPUT.o ^
gfortran_build\obj\Main\SOLVER_MANAGER.o ^
gfortran_build\obj\Main\SUMMARY.o ^
gfortran_build\obj\Main\Initiation\CHECK_DATA.o ^
gfortran_build\obj\Main\Initiation\INIT_ELE_LIB.o ^
gfortran_build\obj\Main\Initiation\INIT_INTEGRATION.o ^
gfortran_build\obj\Main\Initiation\INIT_MECHANICAL.o ^
gfortran_build\obj\Main\Initiation\INIT_SHAPE.o ^
gfortran_build\obj\Main\Initiation\INIT_SOLUTION.o ^
gfortran_build\obj\Main\Initiation\INITIATE.o ^
gfortran_build\obj\Main\Initiation\LOAD_SHAPE.o ^
gfortran_build\obj\Main\Initiation\SHAPE_2D.o ^
gfortran_build\obj\Main\Read\READ_ELEMENTS.o ^
gfortran_build\obj\Main\Read\READ_GEOMETRIES.o ^
gfortran_build\obj\Main\Read\READ_INSTRUCTION.o ^
gfortran_build\obj\Main\Read\READ_MATERIALS.o ^
gfortran_build\obj\Main\Read\READ_XYZ.o ^
gfortran_build\obj\Boundary\BOUNDARY_CONDITION.o ^
gfortran_build\obj\Boundary\BOUNDARY_DISPLACEMENT.o ^
gfortran_build\obj\Boundary\BOUNDARY_FORCE.o ^
gfortran_build\obj\Boundary\READ_DISPLACEMENT.o ^
gfortran_build\obj\Boundary\READ_LINE_FORCE.o ^
gfortran_build\obj\Boundary\READ_NODAL_FORCE.o ^
gfortran_build\obj\Force\GET_FORCE.o ^
gfortran_build\obj\Force\LINE_FORCE_PLANE.o ^
gfortran_build\obj\Mechanics\Stiffness\MATERIAL_2D.o ^
gfortran_build\obj\Mechanics\Stiffness\STIFFNESS.o ^
gfortran_build\obj\Mechanics\Stiffness\STIFFNESS_PLANE.o ^
gfortran_build\obj\Mechanics\Strain\GET_STRAIN.o ^
gfortran_build\obj\Mechanics\Strain\STRAIN_GAUSS_PLANE.o ^
gfortran_build\obj\Mechanics\Strain\STRAIN_NODAL_RECTANGLE.o ^
gfortran_build\obj\Mechanics\Strain\STRAIN_NODAL_TRIANGLE.o ^
gfortran_build\obj\Mechanics\Stress\GET_STRESS.o ^
gfortran_build\obj\Mechanics\Stress\STRESS_GAUSS_PLANE.o ^
gfortran_build\obj\Mechanics\Stress\STRESS_NODAL_RECTANGLE.o ^
gfortran_build\obj\Mechanics\Stress\STRESS_NODAL_TRIANGLE.o ^
gfortran_build\obj\Solver\ASSEMBLE_STIFFNESS.o ^
gfortran_build\obj\Solver\CHOLESCKEY.o ^
gfortran_build\obj\Solver\INV.o ^
gfortran_build\obj\Solver\SOLVE.o ^
gfortran_build\obj\Solver\SOLVER_MECH_STATIC_DIR.o ^
gfortran_build\obj\Main\FEMT.o

echo.
echo Compilation completed. The executable is located at: gfortran_build\FEMT.exe
echo To run the program, use: gfortran_build\FEMT.exe

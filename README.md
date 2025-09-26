# ReDimLungs: Voxel Phantom Organ Matching with Individual Patient Anatomy

This Fortran program semi-automatically modifies a 3D voxel phantom to match the anatomy of an individual patient, specifically by scaling the organs. The script was developed as part of a master thesis on personalized patient dosimetry using Monte Carlo methods and aims to improve the accuracy of organ dose estimates in computational dosimetry, namely Computed Tomography (CT) examinations.

## Description

The ReDimLungs program reads a standard voxel phantom and a user-defined task list to perform volumetric modifications on a specified organ. reDimLungs was developed to work specifically on the lungs. Its primary function is to adjust the dimensions of the left or right lung to match patient-specific measurements. The program's workflow involves:

- Reading configuration parameters and phantom dimensions from an input file: reDimLungs.in.
- Identifying the target organ (right or left lung) and its initial dimensions within the phantom.
- Applying a series of modifications (adding or deleting voxel layers) based on the instructions provided in the tasks file: tasks.dat.
- Writing the modified phantom data to a new output files suitable for visualization with gnuplot and Monte Carlo simulations: ct-den-matXY.dat and phantom.vox.

The work performed with the program employed the Adult Female Reference voxel phantom of the International Commission of Radiological Protection (ICRP-AF), which is a standard model for radiological dose calculations. The modification of the phantom size and lungs can significantly influence the estimated radiation dose imparted to an individual during a CT examination.

## Prerequisites

To compile and run this program, you will need a Fortran compiler. The code has been tested with and is compatible with gfortran.

    Fortran Compiler: A compiler such as gfortran.

    Input Files:

        reDimLungs.in: The primary configuration file.

tasks.dat: A text file containing the list of modifications to apply.

organlist1.dat: A list of organ IDs, material IDs, and densities.

ct-den-matXY_ICRP-AF.dat: The input voxel phantom data file. This file name is defined in 

reDimLungs.in.

## Usage

    Compile the program: Open your terminal or command prompt and compile the source code using gfortran.
    Bash

gfortran -o reDimLungs reDimLungs.f95

This will create an executable file named reDimLungs.

Configure input files: Before running the program, ensure that all the prerequisite input files (reDimLungs.in, tasks.dat, organlist1.dat, and your phantom file) are in the same directory as the executable.

Run the program: Execute the compiled program from the command line.
Bash

./reDimLungs

The program will read the input file and a log of its operations will be written to 

reDimLungs.out.

## Configuration

All program parameters are controlled through the reDimLungs.in file, which uses Fortran's namelist format. The file is separated into distinct sections:

    &PHANTOM_CONFIG: This is the main namelist that encapsulates all configurable parameters.

Phantom & File Settings:

    organ_choice: 1 for the right lung, 2 for the left lung.

phantom_file: The name of the input phantom file.

tasks_file: The name of the file with modification tasks.

tasks_file_lines: The number of task lines to read from the tasks file.

organlist_file: The file containing the list of organs and materials.

organlist_file_lines: The number of lines in the organ list file.

Phantom Dimensions:

    n_slices, n_rows, n_cols: The number of voxels along the z, y, and x axes, respectively.

x_res, y_res, z_res: The resolution (in cm) of each voxel along the x, y, and z axes.

Organ IDs:

    A list of integer IDs for various organs, which must match the IDs in the organlist1.dat file. These IDs are crucial for the program to correctly identify and modify specific organs.

tasks.dat Format

The tasks.dat file specifies the modifications to be applied. Each line represents a single task in the format: 

[axis] [action] [side].

Parameter	Values	Description
axis	x, y, z	The axis to modify
action	a, d	a for add voxel layers, d for delete
side	p, n	p for positive direction, n for negative

Example:

x a p  ! Adds a layer of voxels in the positive x direction
z d n  ! Deletes a layer of voxels from the negative z direction

## Example

The following example shows a typical configuration and the corresponding output for modifying the right lung of the ICRP-AF phantom.

reDimLungs.in

Ini, TOML

&PHANTOM_CONFIG
organ_choice = 1                                 ! 1: Right Lung, 2: Left Lung
phantom_file = 'ct-den-matXY_ICRP-AF.dat'
tasks_file = 'tasks.dat'
tasks_file_lines = 4
organlist_file = 'organlist1.dat'
organlist_file_lines = 142

n_slices = 348         ! Number of slices (z)
n_rows = 137           ! Number of rows (y)
n_cols = 299           ! Number of columns (x)
x_res = 0.1775         ! Voxel resolution in x (cm)
y_res = 0.1775         ! Voxel resolution in y (cm)
z_res = 0.484          ! Voxel resolution in z (cm)

right_lung_id = 45              ! Right Lung Organ ID
...
&END

reDimLungs.out (snippet)

>>>>>>>>>>       ReDimLungs Fortran Routine       >>>>>>>>>> 
...
>>>>> Input file 'reDimLungs.in' read successfully.
The organ to redimension is: Right Lung
...
>>>>> Listing the initial coordinates on all 3 axis
> X AXIS LUNG COORDINATES:
Min x : 54    Max x:245
Size:192 slices; 34.0800 cm.
> Y AXIS LUNG COORDINATES:
Min y : 41    Max y:128
Size: 88 slices; 15.6200 cm.
> Z AXIS LUNG COORDINATES:
Min z :272    Max z :303
Size: 32 slices; 15.4880 cm.
...
>>>>> SECTION: Final Organ Dimensions Analysis.
Recalculate dimensions after all modifications have been applied.
> X AXIS LUNG COORDINATES:
Min x : 41    Max x:240
Size:200 slices; 35.5000 cm.

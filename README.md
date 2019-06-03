# AD2CP_processing
This repository contains step by step instructions on how to process Slocum glider mounted AD2CP data.

This process requires both the dgroup data from the Slocum glider and the AD2CP data from the Nortek AD2CP.

Dependencies:
Slocum Power Tools: https://github.com/kerfoot/spt


The first step in the process is AD2CP_mag_correct.m which corrects the heading data stored in the AD2CP file. This correction is necessary because the AD2CP does not account for the changes in the magnetic field as the pitch battery of the glider moves.

The second step is the AD2CP_beam2enu.m which converts the beam velocities recorded by the AD2CP to East/North/Up direction velocities. There are also some QA/QC steps taken in this program and those thresholds can be adjusted by need.

The third and final step is AD2CP_ls_inversion.m which splits the ENU velocities into components of ocean and glider velocity. This is done by calling a function called inversion_leastSquare_sparse_2019.m which uses a method derived from Visbeck 2002 and Todd et al 2017 to turn a profile of velocity measurements into a system of equations and solve it using a least square inversion. After the inversion, the velocities are split into ocean and glider components and the data is stored in a mat file.

The other files in this repository are tools either required for the processing done in the three steps described above or further data analysis that I have done.

GitHub does not allow storage of large files (understandbly) and these data files can be sizeable. If you are interested in using this package and would like some sample data, please feel free to contact me at sjc244@marine.rutgers.edu and I would be happy to share a sample data set. 

# Pyhsarum-Polycephalum

This repository contains all the scripts used to develop my Bachelor's thesis: "Physarum polycephalum inspired algorithm for the artificial generation of Madrid's metro and train network".  The initial code, and core of this thesis, was obtained from GitHub \url{https://github.com/ammitra/Physarum/blob/main/main.py}. An explanation of each script functionality is given hereafter.
- PPA.py: the main script. It is a OPP algorithm based on Physarum polycephalum. It contains two classes: one to define the Environment and other define the Particles. Each class has several methods. An explanation of the methods' functionality is within the code. It can be compared with the original script. The modifications and effect of the new methods and parameter values can be found in the BEP. 
- pp_image_process.py: script to get the population percentage, pp, corresponding to the metro and train cases. It is based on a RGB filter. This script was run prior to the PPA to determine the pp parameter value. Some methods similar to the code of the script are used within the PPA.py script to localize the FSs and illuminated areas. 
- adjacency_matrix_coefficient.R: script to analyse the adjacency matrices constructed with Microsoft Excel and converted to TSV. The script runs trough the matrices and gets the adjacency matrices values. The script returns Jaccard Index values. For this, it uses some functions from the next script
- la.functions.R: script with several functions for statistical and machine learning analysis. 


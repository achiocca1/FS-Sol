ALL THE FILES PRESENT IN THE FOLDER REFER TO A READY-TO-USE MODEL.

JUST TRY TO RUN FSSol_v1.m IN MATLAB ENVIRONMENT TO TEST THE CODE.

************************
# DESCRIPTION OF THE FILES 
************************

- **ANSYS_Post_Process_Results.mac** = APDL command necessary for a correct post-processing of results in Ansys (it can be used only in Ansys, if any other FEM software is used you need to correctly post-process results as the ones contained in the "FEM_results" folder)

- **FSSol_v1.m** = main Matlab function to be run

- **Results** = folder containing all the nodal results of stress and strain from the finite element analysis

*********************************
# NECESSARY STEPS TO RUN THE SCRIPT 
*********************************

1) **Open FSSol_v1.m**
2) **Set in PARAMETERS the quantities you are interested in:**

- directoryRESULTS      % Directory containing stress and strain results for each node
- blocklength           % How many load steps are present in the .csv files inside "directoryRESULTS"
- ref_times             % Considered loadsteps
- kFS                   % Material constant of Fatemi-Socie critical plane factor
- Sy                    % Yield strengh

3) **RUN the script**
4) **ENJOY!**

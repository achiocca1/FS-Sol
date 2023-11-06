******************************************
DESCRIPTION OF WHAT YOU FIND IN THE FOLDER
******************************************

Results = folder containing all the nodal results of stress and strain from the finite element analysis
FSSol_v1.m = main Matlab function to be run

1) Open FSSol_v1.m
2) Set in PARAMETERS the quantities you are interested in:

- directoryRESULTS      % Directory containing stress and strain results for each node
- blocklength           % How many load steps are present in the .csv files inside "directoryRESULTS"
- ref_times             % Considered loadsteps
- kFS                   % Material constant of Fatemi-Socie critical plane factor
- Sy                    % Yield strengh

3) RUN the script
4) ENJOY!
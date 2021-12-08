# ImprovedEfficiency_via_PlatformDeisgns

Simulation codes for the manuscript entitled "Improved Efficiency for Cross-Arm Comparisons via Platform Designs".

1. The folder *figures* contains all the codes and results (saved in *Rdata* files) that generates all the figures in the manuscript.

2. The file *code_KMbased_Coxbased_InfluFuncs.R* contains the main functions that are used in this simulation, while the file *preliminary_functions.R* collects those needed but basic functions.

3. The file *generate_survivaldata.R* contains the data-generating functions for this simulation. 

4. The file *sim_efficiencygain_of_platformtrial.R* provides the code of a single Monte Carlo run for giving the inference of a given contrast of (unadjusted, stratified and adjusted) relative risks for an arbitrary pair of active interventions relative to the control, which can be used for comparing the efficiency of running platform trials over running multiple separate trials as shown in Section 4 of the manuscript, under the full setting in our submitted manuscript or a selected setting for a quick coding check.

5. The file *sim_adaptive_noninferiority_test_of_platformtrial.R* provides the code of a single Monte Carlo run for giving the rejection rates of the adpative noninferiority test and the likelihood ratio test as shown in Section 5 of the manuscript, under the full setting in our submitted manuscript or a selected setting for a quick coding check.

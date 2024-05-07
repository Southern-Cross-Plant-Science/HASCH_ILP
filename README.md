# HASCH_ILP
SNP panel design using Integer Linear Programming with example application for Cannabis HASCH genotyping platform


This script takes an input conf file and generates and LP file for input to Gurubi or SCIP optimizer.
It can also process the optimizer output for visualization and generation of selected probes list

Usage:
python input.conf soln

Arguments:
input.conf 	 required config file
soln 	 [sol, log] (optional) result file from Gurubi(sol) or SCIP(log) 


The conf file are Python variable assignments, variable=value. Commments should start with #
The variables defined in conf set the execution inputs and optimization options.
Just copy the provided conf and change the values. 
For optional variable, set the value to None, do not remove the assignment

Optional input file flankcountfile containing counts of flanking SNPs for all markers based on strictly filetered vcf. The file has four columns [contig	bp	within100_adj10	within100], 
where within100_adj10 are number of neighbor SNPs within 10bp, and within100 are withn 100bp on both sides. To ignore flanking SNPs, set the filenames to None


The matrix P and pair-wise constraints are saved into files  (*.npy, *.lp.pconstraints) since they are time consuming to generate. Unless recalcP is set to True, the npy and lp.pconstraints files are used on future reruns for the same conf input file. 
Reruns are performed to generate new LP files for varying objectives or weight adjustments. On each settings/objective/constraints/weights used, also change the  'note' variable  below to uniquely identify the output LP file
Reruns are also done to process the optimizer result file. When processing result, the same conf file should be used to match with that of the LP file used.

The optimizer result file should be named  LP_filename + '.sol' for Gurubi, or LP_filename + '.log' for SCIP

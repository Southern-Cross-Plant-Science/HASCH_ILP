# SNP genotyping panel design by Integer Linear Programming
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11149359.svg)](https://doi.org/10.5281/zenodo.11149359)

Used in ***HASCH - A High-throughput Amplicon-based SNP-platform for medicinal Canna-bis and industrial Hemp genotyping applications***

This script takes an input vcf file and generates an LP file for Gurobi or SCIP optimizer. The optimization selects the subset of SNPs that maximizes the number of homozygous mismatches between all sample-pairs in the input genotype data.
It can also process the optimizer output for visualization and generation of selected probes list

**Usage:**
``python vcf2ilp.py input.conf soln``

**Arguments:**

- input.conf  :  required config file
- soln  : (optional) result file from Gurobi(sol) or SCIP(log)  [sol, log]


The conf file are Python variable assignments, variable=value. Commments should start with #
The variables defined in conf set the execution inputs and optimization options.
Just copy the provided conf file and change the values for other inputs. For the optional variables, set the value to ``None`` if not needed.

The optional input file flankcountfile containing counts of flanking SNPs for all markers based on a larger union/semi-filtered SNPs. The file has four columns [contig	bp	within100_adj10	within100], where within100 are number of neighbor SNPs within 100bp, while within100_adj10 are neighbors within 100bp that have neighbors within 10bp. This file can be generated using the included utility scripts countflankingsnps.py/countflankingsnps-bycontig.py. To ignore flanking SNPs, set the flankcountfile to None

The matrix P and pair-wise constraints are saved into files  (*.npy, *.lp.pconstraints) since they are time consuming to generate. Unless recalcP is set to True, the npy and lp.pconstraints files are used on future reruns for the same conf input file. Reruns are performed to generate new LP files for varying objectives or weight adjustments. Each set of settings/objective/constraints/weights should use a unique input.conf, the conf filename is used to uniquely identify the output LP file. Reruns are also performed to process the optimizer result file. When processing result, the same conf file should be used to match with that of the LP file used.


Using the included hasch.conf input, vcf2ilp.py hasch.conf will generate ~75Gb of files (uncompressed) and takes at least 8hrs to run. The generated LP file is accepted by [Gurobi](https://www.gurobi.com/solutions/gurobi-optimizer) or [SCIP](https://www.scipopt.org). The optimizer result file should be named  ``LP_filename + .sol`` for Gurobi, or ``LP_filename + .log`` for SCIP. Using the result file, generate the final probe list and histograms using  ``vcf2ilp.py hasch.conf sol`` or  ``vcf2ilp.py hasch.conf log``




**Required python libraries**
- numpy
- matplotlib
- scipy
- pympler


**Required software**
- plink

 
**Sample inputs files (available in Zenodo repository)**

hasch_input_filtered_57k.vcf.gz

pos57k_ref5m-all.txt


**Utilities**

Generates the flanking SNPs count   (like pos57k_ref5m-all.txt)

``countflankingsnps-bycontig.py``

``countflankingsnps.py``

**Required for countflankingsnps**
- sqlite3
- pandas



**Utility inputs files (available in Zenodo repository)**

cs10-gbs-21trich-wgs7ds-allsnps.merge.maf2_fmis6.vcf.gz-poslist.txt   :  SNP list to count neighborhood (from union of WGS7DS, 21TRICH, GBS then filtered)

hasch_input_filtered_intersect.map  : (snp list for hasch_input_filtered_57k.vcf.gz)




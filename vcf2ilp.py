#
# The MIT License (MIT)
# Copyright © 2024 Locedie Mansueto
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), 
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

'''
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

'''


import sys
import numpy as np
import math
import matplotlib.pyplot as plt

import scipy 
import os
import scipy.sparse as ss
import time

import numpy.matlib as ml
from pympler import asizeof

from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import copy


soln='' 
conffile=sys.argv[1]  # config file
if len(sys.argv)>2:
	soln=sys.argv[2]  # blank:generate LP file, log:read SCIP result, sol:read Gurobi result

with open(conffile) as fin:
	line=fin.readline()
	while line:
		line=line.strip()
		if line and not line.startswith("#"):
			print(line)
			exec(line,globals(),globals())
		line=fin.readline()


#
# INTERNAL FILE NAMING SETTINGS  (just keep these, used in pipeline)
#
fname=vcffile+".plinkA"
contign='allchr'
# used in pipelinings
#contign=sys.argv[4]
#contign=''
if contign:
	fname=fname+'.'+contign

#
# OUTPUT FILES GENERATION SETTINGS
#
#  constraints file
p_constraint_file=fname+'.lp.pconstraints'

# use the conf filename as note,  note is appended to the output LP filename (to document settings defined above) 
note=conffile.replace(".conf","")

# output LP file
foutname=fname+'.' + note + '.lp' 

# flag tasks to do based on command arguments

solnfile=''
if soln=='sol':
	solnfile=foutname + '.sol'
	recalcP=False

elif soln=='log':
	solnfile=foutname + '.log'
	recalcP=False
if solnfile:
	enable_evaluate_pairwise_solution=True

def main():

	#recode vcf into 0,1,2 allele matrix
	if not os.path.exists( fname + '.info.txt') or recalcP:
		print('generating ' + fname + '.info.txt ' +  fname + '.geno.txt')
		create_genoinfo(vcffile)



	# get genomic region intervals indices
	[idx_ends,bp_ends,included_pos,included_poslist,ctg_ends]=get_regions_index(fname)

	#print(str(len(idx_ends)) + 'region boundaries')
	#print("index ends\tctg\tbp ends")
	#for i in range(len(idx_ends)):
	#	print(str(idx_ends[i]) + "\t" + ctg_ends[i] + "\t" + str(bp_ends[i]))
	print('included_pos=' + str(len(included_pos))) # + ' , ' + str(included_pos[0:5]))
	print('interval contraints=' + str(len(idx_ends)))


	# get number of flanking SNPs with low, very low, and high confidence (using SNPs filtered using varying threshold)
	mapPos2Flanks=get_pos2flank_maps(flankcountfile)

	# generate marker weights vector
	obj_weight= [np.nan] * len(included_pos)
	obj_weight=get_objweight_byflanks(obj_weight,included_pos,included_poslist,mapPos2Flanks)
	print('obj_weight=' + str(obj_weight[0:5]) + '...')


	if not enable_evaluate_pairwise_solution:
		# plot weight distrib
		bar_hist(obj_weight, 'objective weight distribution', 'objective fxn weight','No of markers',fname+'.' + note + '.objweight-hist-bar.png')

		print('writing obj_weight histogram objweight-hist.png' )
		plt.hist(obj_weight, bins=100) 
		plt.title('objective weight distribution')
		plt.xlabel('objective fxn weight')
		plt.ylabel('No of markers')
		plt.savefig(fname+'.' + note + '.objweight-hist.png',dpi=600)
		plt.clf()		

		if check_objwt_only:
			exit()


	geno=None
	nrows=None
	ncols=None

	if os.path.exists(fname+'.P.npy') and not recalcP:

		print("P file exists, loading from " + fname+'.P.npy')
		with open(fname + '.P.npy', 'rb') as fnpy:
			sum_p_cols=np.load(fnpy)
			p_row=np.load(fnpy)
			p_col=np.load(fnpy)
			p_data=np.load(fnpy)
			listNAsamples=np.load(fnpy)
			listNBsamples=np.load(fnpy)
			size_geno=np.load(fnpy)
		print('loaded ' + fname + '.P.npy')
		print('size_geno=' + str(size_geno))
		nrows=size_geno[0]
		ncols=size_geno[1]

	else:

		# just get number of samples from plink output, to be used for validations below
		nfams=0
		with open(fname +  '.nosex') as fin:
			line=fin.readline()
			while line:
				nfams+=1
				line=fin.readline()
			print('samples =' + str(nfams))


		print("reading " + fname)
		geno=np.array(list(map(list, np.genfromtxt(fname + '.geno.txt', delimiter='\t',  missing_values="NA", filling_values=np.nan, names=True, dtype=np.float16))))
		nrows=len(geno)
		ncols=np.size(geno,1)
		print(type(geno))
		print('geno rows=' + str(np.size(geno,0)) + ' cols=' + str(np.size(geno,1)))
		#print('rows=' + str(nrows) + ' cols=' + str(ncols(geno)))
		if nfams!=nrows:
			print("fam and geno samples did not match " + str(nfams) + "," + str(nrows))
			exit()

		print('generate and save P matrix to .P.npy,  save P constraints to ' + p_constraint_file)


		# Generate and save P matrix, and sample-pair constraints file (to append into final LP file)
		# p_row,p_col,p_data  are vector representation of matrix P, used for sparse matrix storage, only non-zero elements of P are included
		#p_row  row indices
		#p_col 	column indices
		#p_data data, row and column indices in p_row, and p_col
		# NA and NB matrices where generated from the pairwise combinations of geno rows, such that P=abs(NA-NB) and p_ik final = 0 if p_ik=1 and p_i:=2 for any columnn in i
		#listNAsamples   sample index iof NA matix
		#listNBsamples 	 sample index iof NB matix
		[ p_row,p_col,p_data ,listNAsamples,listNBsamples]=generate_P(geno, nrows, p_constraint_file) 

	print('Non-zero P matrix elements='  + str(len(p_data)) +  '  rows=' + str(len(listNAsamples)))


	pconstset=set()

	#from operator import itemgetter

	print('len(listNBsamples), len(listNAsamples):' + str(len(listNBsamples)) + "," + str(len(listNAsamples)))

	if len(listNBsamples) != len(listNAsamples):
		print('len(listNBsamples) != len(listNAsamples):' + str(len(listNBsamples)) + "," + str(len(listNAsamples)))
		exit()


	if enable_evaluate_pairwise_solution:
		# generate plots and lists from solution markers
		npairs=math.floor(size_geno[0]*(size_geno[0]-1)/2)
		P=ss.coo_matrix( (p_data, (p_row,p_col)), shape=(npairs, size_geno[1]), dtype=np.float16 ) #, dtype=np.byte)
		evaluate_solution(obj_weight,included_poslist,sum_p_cols,size_geno,P,listNAsamples,listNBsamples)
		exit()


	# generate LP file
	print('Generating LP file')


	with open(foutname,'w',buffering=-1) as fout:
		xs=[]
		xscount=[]
		xsbound=[]
		xsbinaries=[]

		prows= math.floor(nrows*(nrows-1)/2) 


		# 
		#	generate objective function,
		#   and generate boundary and binary constraints,
		#

		CdotP=None
		if enable_maximize_mismatches:
			# calculate cT P for objective function
			#CdotP = cT x P in objective functon, c is 1 vector
			C = ss.coo_matrix( (  [1]*prows, ([0]*prows, range(prows))), shape=(1,prows) ) 
			P = ss.coo_matrix( (p_data, (p_row,p_col)), shape=( prows,ncols), dtype=np.float16) 
			CdotP=C.tocsr().dot(P.tocsr())	
			print('size CdotP=' + str(np.size(CdotP,0)) + ' x ' + str(np.size(CdotP,1)))

			print('Generating Maximize mismatches objective')
			fout.write("Maximize\n obj: ")
		else:
			print('Generating Minimize markers objective')
			fout.write("Minimize\n obj: ")

		npairs=prows
		nanCdotPcount=0
		zeroCdotPcount=0
		for i in range(1,ncols+1):
			if enable_maximize_mismatches:
				cdotpij=CdotP[0,i-1]
				if np.isnan(cdotpij):
					nanCdotPcount+=1
				elif cdotpij==0.0:
					zeroCdotPcount+=1
				else:
					if flankcountfile:
						# adjust using weight if provided
						xs.append( "{:.3f}".format(obj_weight[i-1]*cdotpij) + ' x' + str(i))
					else:
						xs.append( "{:.3f}".format(cdotpij) + ' x' + str(i))
			else:
				xs.append( 'x' + str(i))

			#  for objective equation
			xscount.append('x' + str(i) )
			
			# for boundary and binaryconstraints
			xsbinaries.append('x' + str(i))
			xsbound.append('0 <= ' + 'x' + str(i) + ' <= 1')

		fout.write(' + '.join(xs) + "\n")

		print('nanCdotPcount=' + str(nanCdotPcount))
		print('zeroCdotPcount=' + str(zeroCdotPcount))

		fout.write("Subject To\n")

		conscount=0
		xs=[]


		# 
		#	generate constraints equations
		#

		global maxmarkers 	#  maxmarkers can change if adjust_maxmarkers is set
		if enable_maximize_mismatches:
			print('Generating maximum markers constraint')
			if len(idx_ends)>maxmarkers:
				print('len(idx_ends)>maxmarkers:' + str(len(idx_ends)) + '  > ' + str(maxmarkers))
				if adjust_maxmarkers:
					print('adjusted maxmarkers to:' + str( math.floor(len(idx_ends)*adjust_maxmarkers)))
					maxmarkers=math.floor(len(idx_ends)*adjust_maxmarkers)
				else:
					exit()
			conscount+=1
			fout.write(" c1: " + " + ".join(xscount) + " <= " + str(maxmarkers) + "\n")

		if use_bp_intervals:
			print('Generating bp regions constraints')
			print('intervals=' + str(len(idx_ends)))
			idx_end_count=0
			next_idx_end=idx_ends[idx_end_count]
			interval_constraints_cnt=0
			for i in range(1,ncols+1):
				xs.append('x' + str(i))
				if i>next_idx_end:
					conscount+=1
					interval_constraints_cnt+=1
					fout.write(' c' + str(conscount) + ': ' + ' + '.join(xs) + " >= " + str(minperinterval) + "\n")

					xs=[]
					idx_end_count+=1
					if idx_end_count<len(idx_ends):
						next_idx_end=idx_ends[idx_end_count]
					else:
						next_idx_end=1000000000
			print('interval_constraints_cnt=' + str(interval_constraints_cnt))
		else:
			print('Generating marker index constraint')
			for i in range(1,ncols+1):
				xs.append('x' + str(i))
				if i % maxinterval==0:
					conscount+=1
					fout.write(' c' + str(conscount) + ': ' + ' + '.join(xs) + " >= " + str(minperinterval) + "\n")
					xs=[]

		if xs:
			conscount+=1
			fout.write(' c' + str(conscount) + ': ' + ' + '.join(xs) + " >= " + str(minperinterval) + "\n")


		print('Appending sample-pair constraints file generated during P calculation ... ' +  p_constraint_file)
		with open(p_constraint_file) as pin:
			line=pin.readline().rstrip()
			pair_constraints_cnt=0
			while line:
				conscount+=1
				pair_constraints_cnt+=1
				fout.write(' c' + str(conscount) + ': ' + line + "\n")
				line=pin.readline().rstrip()
			print('pair_constraints_cnt=' + str(pair_constraints_cnt))
	

		print('bounds contraints=' + str(len(xsbound)))
		print('binary contraints=' + str(len(xsbinaries)))
		fout.write("Bounds\n")
		fout.write('\n'.join(xsbound))
		fout.write("\n")
		fout.write("Binaries\n")
		fout.write(" ".join(xsbinaries))
		fout.write("\nEnd\n")


	print('generated ' + foutname)


def create_genoinfo(vcf):

	# generate map
	#cmd="plink -vcf " + vcf + " --recode -out " + vcf + ".plinkA --allow-extra-chr --double-id --vcf-idspace-to '_'"
	cmd="plink -vcf " + vcf + " --make-bed -out " + fname  +" --allow-extra-chr --double-id --vcf-idspace-to '_'"
	print(cmd)
	os.system(cmd)

	# generate geno, info files
	cmd="plink -vcf " + vcf + " --recode A -out " + fname+ " --allow-extra-chr --double-id  --vcf-idspace-to '_'"
	print(cmd)
	os.system(cmd)
	cmd="cut -d ' ' --output-delimiter '\t' -f 7- " + fname + ".raw | tail -n +2 > " + fname+ ".geno.tmp" 
	print(cmd)
	os.system(cmd)
	cmd="head -n 1 " + fname + ".raw | cut -d ' ' --output-delimiter '\t' -f 7- |  sed 's/^S//g;s/\tS/\t/g;s/_[ATGC]//g' > " + fname + ".header.tmp" 
	print(cmd)
	os.system(cmd)
	cmd="cat " +  fname+ ".header.tmp " +  fname + ".geno.tmp > " + fname + ".geno.txt" 
	print(cmd)  
	os.system(cmd)
	os.system("rm " + fname + ".header.tmp" )
	os.system("rm " +fname + ".geno.tmp" )
	os.system("rm " +fname + ".raw" )

	curChrpos=[]
	with open(fname+ '.bim') as fin, open( fname + '.info.txt',"w") as foutinfo:
		foutinfo.write("chrN\trsId\tbp\n")
		line=fin.readline()
		while line:
			cols=line.rstrip().split("\t")
			curChrpos.append(cols[0]+"_"+cols[3])


			foutinfo.write(cols[0]+"\t" + cols[0]+"_"+cols[3] + "\t" + cols[3] + "\n")
			line=fin.readline()




def create_matrix(n,m,val, scipy=False, dtype=np.float16):
	
	if scipy:
		a=csr_matrix((n, m), dtype= np.ubyte)
		m=csr_matrix((n, m), dtype= np.bool)
		return [a,m]
	else:
		a = np.empty((n,m),dtype=dtype)
		a[:] = val
		return a

def bar_hist(values, title, xlabel,ylabel,filename):
	cnt2val=dict()
	for x in values:
		if x in cnt2val:
			cnt2val[x]+=1
		else:
			cnt2val[x]=1
	distx=[]
	disty=[]
	for xy in list(cnt2val.keys()):
		distx.append(xy)
		disty.append(cnt2val[xy])

	print('writing histogram ' + title + ' ' + filename )
	plt.scatter(distx, disty) 
	#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
	#plt.yscale('log', nonpositive='clip')
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.savefig(filename,dpi=600)
	plt.clf()		


def get_regions_index(fname):

	included_pos=dict()
	included_poslist=[]
	idx_ends=[] # 0-based
	bp_ends=[]
	ctg_ends=[]

	last_idx_end=-1
	prev_bp=None
	prev_ctg=''
	next_bp_end=bp_per_segment

	print('using orig info format')
	with open(fname + ".info.txt") as fin:
		#chrN    rsId    bp
		#NC_029855.1     NC_029855.1_122194      122194
		#NC_029855.1     NC_029855.1_195292      195292
		line=fin.readline()
		line=fin.readline()
		idx_count=0
		while line:
			line=line.rstrip()
			try:
				col2=line.split("\t")
				ctg=col2[0]
				bp=int(col2[2])

				included_pos[col2[1]]=True
				included_poslist.append(col2[1])
				if (ctg==prev_ctg) or not prev_ctg:
					if bp>next_bp_end:
						if last_idx_end<idx_count-1:
							idx_ends.append(idx_count-1)
							bp_ends.append(prev_bp)
							ctg_ends.append(ctg)
							last_idx_end=idx_count-1
						else:
							print('no idx between ' + str(prev_bp) + '-' + str(bp))
						next_bp_end+=bp_per_segment
				else:
					print('ending chromosome ' +  prev_ctg + ' ' + str(prev_bp) + ' idx ' + str(idx_count-1))
					print('starting new chromosome boundaries for ' + ctg + ' ' + str(bp) + ' at idx ' + str(idx_count))
					idx_ends.append(idx_count-1)
					bp_ends.append(prev_bp)
					ctg_ends.append(prev_ctg)
					next_bp_end=bp_per_segment

				idx_count+=1
				prev_bp=bp
				prev_ctg=ctg
			except Exception as e:
				print(str(e))
				print(cols)
				print(col1)
				print(line)
				raise e
			line=fin.readline().rstrip()

		print('ending chromosome ' +  prev_ctg + ' ' + str(prev_bp) + ' idx ' + str(idx_count-1))
		idx_ends.append(idx_count-1)
		bp_ends.append(prev_bp)
		ctg_ends.append(prev_ctg)

	return [idx_ends,bp_ends,included_pos,included_poslist,ctg_ends]



def getSelectedIndex(contig):
	# gete gelected indecces from optimazation soluution
	print('reading solution ' + solnfile)
	includeIdx=[] 
	obj_value=None

	if soln=='log':	# SCIP output
		with open(solnfile) as fin:
			line=fin.readline()
			while not line.startswith('objective value:'):
				line=fin.readline()

			obj_value=line.split(":")[1].strip()
			
			line=fin.readline()
			inx=False
			while line:
				if line.startswith('x'):
					cols=line.split(" ",1)
					idx=int(cols[0].replace("x",""))-1
					includeIdx.append(idx)
					inx=True
				else:
					if inx:
						break
				line=fin.readline()
	elif soln=='sol': 	# Gurubi result
		with open(solnfile) as fin:
			line=fin.readline()

			while not line.startswith('# Objective value'):
				line=fin.readline()
			print(line)
			obj_value=line.split("=")[1].strip()
			line=fin.readline()
			#print(line)
			inx=False
			while line:
				if line.startswith('x'):
					cols=line.split(" ",1)
					if cols[1].strip()=="1":
						idx=int(cols[0].replace("x","").strip())-1
						includeIdx.append(idx)
					inx=True
				else:
					if inx:
						break
				line=fin.readline()
	print('obj_value=' + obj_value + ' non-zero=' + str(len(includeIdx)))
	if not float(obj_value)>0:
		print('Negative objective!')
		#exit()
	return includeIdx



def get_pos2flank_maps(flankcountfile):
	mapPos2Flanks=dict()
	if flankcountfile:
		with open(flankcountfile) as fin:
			line=fin.readline().rstrip()
			#contig	bp	within100_adj10	within100
			line=fin.readline().rstrip()
			i=0
			while line:
				cols=line.split("\t")
				cgtbp=cols[0]+"_"+cols[1]
				mapPos2Flanks[cgtbp]=[float(cols[3]), float(cols[2])]
				line=fin.readline().rstrip()
				i+=1

			print(str(i) + ' flankcount values added')
	return mapPos2Flanks


def get_objweight_byflanks(obj_weight,included_pos,included_poslist,mapPos2Flanks):
	k=0
	for i in range(len(included_pos)):
		ctgbp=included_poslist[i]
		if ctgbp in mapPos2Flanks:
			cntadj=mapPos2Flanks[ctgbp]
			obj_weight[i]=max(maxflanks_value, 1.0-wtflank[0]*cntadj[0]-wtflank[1]*cntadj[1])
			k+=1
		else:
			obj_weight[i]=1

		if obj_weight[i]>1:
			print('obj_weight[i]>1 ' + str(obj_weight[i]) + ' for ' + ctgbp)
			if ctgbp in mapPos2Flanks:
				print('mapPos2Flanks=' + str(mapPos2Flanks[ctgbp]))

	print(str(k) + '/'  + str(len(included_pos)) + ' flankcounts in objective function')
	return obj_weight



def generate_P(geno, nrows, p_constraint_file): 

	hom_pairs_count=0
	het_pairs_count=0
	no_pairs_count=0
	hom_pairs_count_uniq=0
	het_pairs_count_uniq=0
	hom_unchecked_uniq_count=0
	het_unchecked_uniq_count=0

	listNAsamples=[]
	listNBsamples=[]

	p_row=np.array([],dtype=np.uint32)
	p_col=np.array([],dtype=np.uint32)
	p_data=np.array([],dtype=np.uint8)


	geno_cols=np.size(geno,1)
	sum_p_cols=[0] * geno_cols
	Prow2Pairs=[None] * math.floor(nrows*(nrows-1)/2)

	start_time = time.time()
	foutP=open(p_constraint_file,'w',buffering=-1)

	for i in range(nrows):

		# for monitoring progress
		if i % 100==0:
			end_time = time.time()
			print('\nrow=' + str(i) + ":" + str(end_time - start_time))
			start_time=end_time


		istart=len(listNBsamples)
		repi = ml.repmat(geno[i],nrows-i-1,1)
		repirows=istart + (nrows-i-1)

		gcols=np.size(geno,1)
		colstart=0
		colend=colstart+split_col
		colend_actual=min(colend, gcols)

		# repi is repeat of row geno[i], repeated nrows-i-1 times
		# generate nmim=Na-Nb by appending rows of abs(geno[i+1,]-repi[,]),  for each row i in geno G,  NaN (resulting from missing values 'NA') are set to 0
		nmin=np.nan_to_num(np.abs(geno[i+1:,colstart:colend_actual]-repi[:,colstart:colend_actual])).astype(np.byte)

		colstart=colend
		colend=colend+split_col

		# grow columns of nmin split_col columns at a time, (for mmemory mnagement)
		while colstart<=gcols:
			colend_actual=min(colend, gcols)
			nmin=np.hstack( (nmin, np.nan_to_num(np.abs(geno[i+1:,colstart:colend_actual]-repi[:,colstart:colend_actual])).astype(np.byte)))
			colstart=colend
			colend=colend+split_col

		# for monitoring memory usage 
		if i % 100 ==0:
			print("Checking memory usage at i=" + str(i)+ "/" + str(nrows))
			print("size nmin\t" + asizeof.asized(nmin).format())
			print('var dir getsizeof')
			
			
			varsize=[]
			for var in dir():
				#print(str([var, type(eval(var)), eval(var), sys.getsizeof(eval(var)) ] ))
				varsize.append((var, type(eval(var)), sys.getsizeof(eval(var))))
			
			#for var in varsize.sort(key=itemgetter(2)):
			for var in sorted(varsize, key=lambda tup: tup[2], reverse=True):
				print(str(var))


		conslist=[]

		# iterate onky on non-zero NA-NB
		nz_rowcol = np.nonzero(nmin)
		nz_row=nz_rowcol[0]
		nz_col=nz_rowcol[1]
		prevri=0
		mini_hom_const=[]
		mini_het_const=[]
		nz_rows=len(nz_row)
		nz_rows_m1=nz_rows-1


		p_row_i=[] #np.array([], np.uint32)
		p_col_i=[] #np.array([], np.uint32)
		p_data_i=[] #np.array([], np.uint8)
		for nzi in range(nz_rows):
			ri=nz_row[nzi]
			if ri<prevri:
				print('decreasing r in nz_row ' + str(ri) +'<' + str(prevri))
				exit()

			# if new sample-pair ri or last row, generate sample-pair constraint line 
			if ri>prevri or nzi==nz_rows_m1:
				if  nzi==nz_rows_m1:  # if last row, include 
					# read last element before precessing (same as last lines of this for loop)
					ci=nz_col[nzi]
					nminij=nmin[ri,ci]
					p_row_i.append(istart+ri)
					p_col_i.append(ci)
					p_data_i.append( nminij)

					if nminij>1:
						mini_hom_const.append(  '2 x' + str(ci+1))
						sum_p_cols[ci]+=2
					elif nminij==1:
						mini_het_const.append(  'x' + str(ci+1))
						sum_p_cols[ci]+=1


				# new geno row or end, process prev block of nmin=NA-Nb (geno[i+1-end] and repmat[i] rows)
				if len(mini_hom_const)>0:  # pair has homozygous mismatch
					hom_pairs_count+=1
					cons_str=" + ".join(mini_hom_const) + " >= " + str(minperpair) 
					if len(cons_str) > buffer_p_tofile_len:
							foutP.write(cons_str+"\n")
							hom_unchecked_uniq_count+=1
					else:
						if not cons_str in pconstset:  # exclude short, duplicated constraits
							pconstset.add(cons_str)
							hom_pairs_count_uniq+=1
							foutP.write(cons_str+"\n")

				elif len(mini_het_const)>0: # pair has heterozygous mismatch only
					het_pairs_count+=1
					cons_str=" + ".join(mini_het_const) + " >= " + str(minperpair) 
					if len(cons_str) > buffer_p_tofile_len:
							foutP.write(cons_str+"\n")
							het_unchecked_uniq_count+=1
					else:
						if not cons_str in pconstset:  # exclude short, duplicated constraits
							pconstset.add(cons_str)
							het_pairs_count_uniq+=1
							foutP.write(cons_str+"\n")
					hetconstraint_pairs.append(str(i)+','+str(i+ri+1))
				else:
					no_pairs_count+=1
			

				mini_hom_const=[]
				mini_het_const=[]
			
			if nzi<nz_rows_m1:
				ci=nz_col[nzi]
				nminij=nmin[ri,ci]
				p_row_i.append(istart+ri)
				p_col_i.append(ci)
				p_data_i.append(nminij)

				if nminij>1:
					mini_hom_const.append(  '2 x' + str(ci+1))
					sum_p_cols[ci]+=2
				elif nminij==1:
					mini_het_const.append(  'x' + str(ci+1))
					sum_p_cols[ci]+=1
	
			prevri=ri


		p_row= np.append(p_row, np.array(p_row_i, np.uint32))
		p_col= np.append(p_col, np.array(p_col_i, np.uint32))
		p_data= np.append(p_data, np.array(p_data_i, np.uint8))

		for ij in range(istart,repirows):
			Prow2Pairs[ij]=[i, i+1+(ij-istart)]

		listNBsamples=listNBsamples + [i]*(nrows-i-1)
		listNAsamples=listNAsamples + list(range(i+1,nrows))

	foutP.close()

	size_geno=np.shape(geno)

	print('saving P matrix using numpy save')
	try:
		with open(fname + '.P.npy', 'wb') as fnpy:
			np.save(fnpy, np.array(sum_p_cols))
			np.save(fnpy, p_row)
			np.save(fnpy, p_col)
			np.save(fnpy, p_data)
			np.save(fnpy, np.array(listNAsamples))
			np.save(fnpy, np.array(listNBsamples))
			np.save(fnpy, np.array(size_geno))
		print('saved ' + fname + '.P.npy')
	except Exception as e2:
		print(str(e2))
		raise e2

	print('P matix generation statistcs')
	print('geno size ' + str(size_geno))
	print('hom_pairs_count=' + str(hom_pairs_count))
	print('hom_pairs_count_uniq=' + str(hom_pairs_count_uniq))
	print('het_pairs_count=' + str(het_pairs_count))
	print('het_pairs_count_uniq=' + str(het_pairs_count_uniq))
	print('no_pairs_count=' + str(no_pairs_count))
	print('hom_unchecked_uniq_count='+ str(hom_unchecked_uniq_count))
	print('het_unchecked_uniq_count=' + str(het_unchecked_uniq_count))


	return [p_row,p_col,p_data,listNAsamples,listNBsamples]



# generate plots using selected markers
def evaluate_solution(obj_weight,included_poslist,sum_p_cols,size_geno,P,listNAsamples,listNBsamples):

	cols_geno=size_geno[1]
	rows_geno=size_geno[0]
	selIdx=getSelectedIndex(contign)
	obj_weight_sol=[]
	negwt=0
	for ictgbp in selIdx:
		if not obj_weight[ictgbp]>0:
			negwt+=1
		obj_weight_sol.append(obj_weight[ictgbp])
	if negwt>0:
		print(str(negwt) + ' negative objective weights in solution')
		#exit()
	print('selected markers=' + str(len(selIdx)))

	with open(solnfile + '.poslist.txt','w') as solposlist, open(solnfile + '.objwt_sumpcol_list.txt','w') as solobjwtlist:
		for ictgbp in selIdx:
			col1=included_poslist[ictgbp]
			col1=col1.split("1_")
			ctg=col1[0]+"1"
			solposlist.write(ctg + "\t" + col1[1]+"\n")
			solobjwtlist.write(ctg + "\t" + col1[1]+"\t" + str(obj_weight[ictgbp]) + "\t" + str(sum_p_cols[ictgbp]) + "\n")


	print('writing obj_weight histogram of solution solm-objweight-hist.png' )
	plt.hist(obj_weight_sol, bins=100) 
	#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
	#plt.yscale('log', nonpositive='clip')
	plt.title('solution objective weight distribution')
	plt.xlabel('objective fxn weight')
	plt.ylabel('No of markers')
	plt.savefig(solnfile + '.objweight-hist.png',dpi=600)
	plt.clf()		


	if poslist_only:
		exit()


	nfams=0
	with open(fname +  '.nosex') as fin:
		line=fin.readline()
		while line:
			nfams+=1
			line=fin.readline()
		print('samples =' + str(nfams))

	Mismatches=[]
	Nmismatches=create_matrix(nfams,nfams,0)
	MismatchesPerMarker=[]
	zeropairs=[]
	listpairs=[]
	lowestcount=100000000000
	highestcount=0
	paircount_by_mismatch=dict()


	idx_row=[]
	idx_data=[]
	idx_col=[0]* len(selIdx)
	for idx in selIdx:
		idx_row.append(idx)
		idx_data.append(1)

		MismatchesPerMarker.append(sum_p_cols[idx])

	M = ss.coo_matrix( (idx_data, (idx_row,idx_col)), shape=(cols_geno,1)) #, dtype=np.byte) 

	


	# PdotM is number of mismatches for each sample pair at the selected markers M
	PdotM= P.tocsr().dot(M.tocsr())

	print('PdotM=' + str(np.size(PdotM,0)) + ' x ' +  str(np.size(PdotM,1)) ) 
	print(PdotM)
	#print('listNAsamples' + str(listNAsamples))
	#print('listNBsamples' + str(listNBsamples))

	mismatchpairs=set()
	for ij in range(PdotM.shape[0]):
		Pval=PdotM[ij,0]

		if Pval>0:				
			Nmismatches[listNAsamples[ij],listNBsamples[ij]]=Pval
			Nmismatches[listNBsamples[ij],listNAsamples[ij]]=Pval
			if Pval<lowestcount:
				lowestcount=Pval
			if Pval>highestcount:
				highestcount=Pval

			if Pval in paircount_by_mismatch:
				paircount_by_mismatch[Pval]+=1
			else:
				paircount_by_mismatch[Pval]=1

			Mismatches.append(Pval)
			#ijpair=str(listNAsamples[ij]) + "\t" + str(listNBsamples[ij])
			ijpair=str(listNBsamples[ij]) + "\t" + str(listNAsamples[ij])
			listpairs.append(ijpair + "\t" + str(Pval))
			mismatchpairs.add(ijpair)

	for ij in range(size_geno[0]):
		for jk in range(ij+1, size_geno[0]):
			ijpair=str(ij) + "\t" + str(jk)
			if ijpair in mismatchpairs:
				pass
			else:
				zeropairs.append(ijpair)


	print('lowest mismatch count=' + str(lowestcount))
	print('highest mismatch count=' + str(highestcount))
	print('zeros=' + str(len(zeropairs)))

	if len(zeropairs)>0:
		print('writing zero mismatch pairs file ' + solnfile + '.nomatched-bypairs.txt')
		with open(solnfile + '.nomatched-bypairs.txt', 'w') as fout:
			fout.write("\n".join(zeropairs))


	print('writing counts list file ' + solnfile + '.mismatchcounts-bypairs.txt')
	with open(solnfile + '.mismatchcounts-bypairs.txt', 'w') as fout:
		fout.write("\n".join(listpairs))

	#print('writing counts list file ' + solnfile + '.mismatchcounts.txt')
	#with open(solnfile + '.mismatchcounts.txt', 'w') as fout:
	#	fout.write("\n".join(PdotM.tolist()))
	
	#from fast_histogram import histogram1d, histogram2d
	#histogram1d(x, y, range=[[-1, 2], [-2, 4]], bins=30)

	paircount_by_mismatch_value=[]
	for k in paircount_by_mismatch:
		paircount_by_mismatch_value.append(paircount_by_mismatch[k])

	if False:  # plotting

		if False: # reuqires to pass p_row,p_col, p_data in parameter
			print('writing counts histogram ' + solnfile + '.p_col.png')
			plt.hist(p_col, bins=100) 
			#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
			plt.title('P col distribution')
			plt.xlabel('p_col idx')
			plt.ylabel('No entries')
			plt.savefig(solnfile + '.p_col.png',dpi=600)
			plt.clf()		

			print('writing counts histogram ' + solnfile + '.p_row.png')
			plt.hist(p_row, bins=100) 
			#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
			plt.yscale('log', nonpositive='clip')
			plt.title('P row distribution')
			plt.xlabel('p_row idx')
			plt.ylabel('No entries')
			plt.savefig(solnfile + '.p_row.png',dpi=600)
			plt.clf()		

			print('writing counts histogram ' + solnfile + '.p_data.png')
			plt.hist(p_data, bins=100) 
			#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
			plt.title('P data distribution')
			plt.xlabel('p data')
			plt.ylabel('No entries')
			plt.savefig(solnfile + '.p_data.png',dpi=600)
			plt.clf()		

		print('writing counts histogram ' + solnfile + '.mismatch_by_pair_bar.png')
		plt.bar(list(paircount_by_mismatch.keys()), paircount_by_mismatch_value) 
		#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
		plt.title('Mismatch counts distribution (bar)\n0s=' + str(len(zeropairs)))
		plt.xlabel('No of allele mismatches')
		plt.ylabel('No of sample pairs')
		plt.savefig(solnfile + '.mismatchbypair-histogram_bar.png',dpi=600)
		plt.clf()		

		if True:

			print('writing counts histogram ' + solnfile + '.mismatch_by_pair_histogram.png')
			plt.hist(Mismatches, bins=100) 
			#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
			plt.title('Mismatch counts distribution (hist)\n0s=' + str(len(zeropairs)) + ' markers=' + str(len(selIdx)))
			plt.xlabel('No of allele mismatches')
			plt.ylabel('No of sample pairs')
			plt.savefig(solnfile + '.mismatchbypair-histogram.png',dpi=600)
			plt.clf()		

			print('writing counts histogram ' + solnfile + '.mismatch_by_pair_histogram-log.png')
			plt.hist(Mismatches, bins=100) 
			#plt.hist(Mismatches, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, *, data=None, **kwargs)
			plt.yscale('log')
			plt.title('Mismatch counts distribution\n0s=' + str(len(zeropairs))+ ' markers=' + str(len(selIdx)))
			plt.xlabel('No of allele mismatches')
			plt.ylabel('No of sample pairs')
			plt.savefig(solnfile + '.mismatchbypair-histogram-log.png',dpi=600)
			plt.clf()		

		if True:
			print('writing count mismatch per marker ' + solnfile + '.mismatch_by_marker.png')
			#fig, ax = plt.subplots()
			plt.bar(range(size_geno[1]), sum_p_cols, color='cyan', alpha=0.5,linewidth=0) 
			plt.bar(selIdx, MismatchesPerMarker,color='red',width=1.2,linewidth=0) 
			plt.xlim(1, size_geno[1])
			plt.xlabel("Marker no.")
			plt.ylabel("No. of allele mismatches")
			plt.savefig(solnfile + '.mismatchbymarker-bar.png',dpi=600)


			print('writing count mismatch per marker ' + solnfile + '.mismatch_by_marker-log.png')
			#fig, ax = plt.subplots()
			plt.bar(range(size_geno[1]), sum_p_cols, color='cyan', alpha=0.5,linewidth=0) 
			plt.bar(selIdx, MismatchesPerMarker,color='red',width=1.2,linewidth=0) 
			plt.yscale('log')
			plt.xlim(1, size_geno[1])
			plt.xlabel("Marker no.")
			plt.ylabel("No. of allele mismatches")
			plt.savefig(solnfile + '.mismatchbymarker-bar-log.png',dpi=600)

	exit()


main()
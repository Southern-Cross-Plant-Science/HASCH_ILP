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

# This script counts the number of neighbors for the SNPs in poslist, using the neighborhood/reference in reflist
# First edit the input files, labels in the variables below

# USAGE:
# python countflanking.py
# python countflanking.py include_chr
# python countflanking.py "" exclude_chr

# The ouput file is posname_refname-chr.txt or posname_refname-all.txt, 
# posname,refname are defined below

# The output is table with 4 columns [contig,bp,within100_adj10,within100]
# within100 is the number of neighbors within 100bp on both sides
# within100_adj10 is the number of neighbors within 100bp on both sides, that have neighbors within 10bp
#
# If run by chromosome, manually concatenate all outputs to get the final results 
# This script can also called from countflankingsnps-bycontig.py to automatically start job per chromosome
#


import sys
import os
import sqlite3 as sl
import pandas as pd
import csv


# run only this chromosomes with name containing
includechr=''
# exclude contigs with name containing
excludechr=''

# optionally set above from command argument
if len(sys.argv)>1:
	includechr=sys.argv[1]
if len(sys.argv)>2:
	excludechr=sys.argv[2]


# list with chromosome and bp columns of SNPs to count neighbors
poslist='hasch_input_filtered_intersect.map'
# column number of chromosome and bp (1-index)
poscols='1,4'
# name of snplist
posname='pos57k'

# list with chromosome and bp columns where neighborhood is counted

#reflist='cs10-ipkgbs-21trichs-wgs7ds-allsnps.merge.filltags.maf2_fmis6.vcf.gz-poslist.txt'
reflist='cs10-gbs-21trich-wgs7ds-allsnps.merge.maf2_fmis6.vcf.gz-poslist.txt' 
# name of neighbors snplist
refname='ref5m'
# column number of chromosome and bp (1-index). if None, assumed only two columns at 1,2
refcols=None


# count distance pairs by batch of this bp length
byrange=10000
# regenerate database
recalc=True
# testing only to show SQLs
printsqlonly=False



dbname=posname + '_' + refname + '.db'
if includechr and excludechr:
	print('set includechr or excludechr only, not both')
	exit()
if includechr:
	dbname=posname + '_' + refname + '-' + includechr + '.db'
	cmd='grep ' + includechr + ' ' +  reflist + ' > ' + reflist + '-' + includechr
	print(cmd)
	os.system(cmd)
	reflist=reflist + '-' + includechr
	cmd='grep ' + includechr + ' '  + poslist + ' > ' + poslist + '-' + includechr
	print(cmd)
	os.system(cmd)
	poslist=poslist + '-' + includechr
if excludechr:
	dbname=posname + '_' + refname + '-' + excludechr + '.db'
	cmd='grep -v ' + excludechr + ' ' +  reflist + ' > ' + reflist + '-' + excludechr
	print(cmd)
	os.system(cmd)
	reflist=reflist + '-' + excludechr
	cmd='grep -v ' + excludechr + ' '  + poslist + ' > ' + poslist + '-' + excludechr
	print(cmd)
	os.system(cmd)
	poslist=poslist + '-' + excludechr


if os.path.exists(dbname):
	print(dbname + " exists..")
	if recalc:
		print("recalculating")
		os.system('rm ' + dbname)
	else:
		print("keeping")


con = sl.connect(dbname)
cur = con.cursor()


def create_table(tablename,fin,inrefcols=None):
	# create table from fin file with chr,bp column numbers in inrefcols
	os.system('rm ' +  fin + '.tmp')
	if inrefcols:
		os.system('cut -f' + inrefcols + ' ' + fin + ' > ' + fin + '.tmp')
	else:
		os.system('ln -s ' + fin +  ' ' + fin + '.tmp')

	cur.execute('CREATE TABLE IF NOT EXISTS ' + tablename + ' (id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, contig TEXT, bp INTEGER);')

	print('loading ' + fin + '.tmp')
	os.system('head -n 10 ' + fin + '.tmp')
	print()
	with open(fin + '.tmp') as csv_file:
		csv_reader = csv.reader(csv_file, delimiter='\t')
		cur.executemany("INSERT INTO " + tablename + "(CONTIG,BP) VALUES(?,?)", csv_reader)
		#con.commit()
	cur.execute("CREATE INDEX IF NOT EXISTS index_" + tablename + "_contigbp ON " + tablename + "(contig, bp);")
	cur.execute("CREATE INDEX IF NOT EXISTS index_" + tablename + "_contig ON " + tablename + "(contig);")
	cur.execute("CREATE INDEX IF NOT EXISTS index_" + tablename + "_bp ON " + tablename + "(bp);")
	print('loaded ' + fin + ' to ' +  tablename + ' in file ' + dbname)


def create_counts_table(postable,reftable,limit):
	# count neigbors between postable and reftable within limit bp
	limit=str(limit)
	cur.execute('CREATE TABLE IF NOT EXISTS ' + postable + '_' + reftable + '_count' + limit +  ' as select pk.contig, pk.bp, count(rk.bp-pk.bp) as count_' + limit + ' from ' + postable + ' pk, ' + reftable + ' rk where pk.contig=rk.contig and abs(rk.bp-pk.bp)<=' + limit + ' and rk.bp<>pk.bp group by pk.contig,pk.bp;')
	cur.execute("CREATE INDEX IF NOT EXISTS index_" + postable + '_' + reftable + '_count' + limit   + "_contigbp ON " + postable + '_' + reftable + '_count' + limit  + "(contig, bp);")
	cur.execute("CREATE INDEX IF NOT EXISTS index_" + postable + '_' + reftable + '_count' + limit   + "_contig ON " + postable + '_' + reftable + '_count' + limit  + "(contig);")
	cur.execute("CREATE INDEX IF NOT EXISTS index_" + postable + '_' + reftable + '_count' + limit   + "_bp ON " + postable + '_' + reftable + '_count' + limit  + "(bp);")
	print('loaded ' + postable + '_' + reftable + '_count' + limit)	

def create_counts_table_byrange(postable,reftable,limit,byrange=None,withintable=None,withinlimit=None):
	# count neigbors between postable and reftable within limit bp
	# if withintable is provided, count the number of SNPs in withintable within withinlimit, where postable and reftable is within limit bp

	limit=str(limit)
	counttable=None

	if byrange:
		sql='select contig, max(bp) as maxbp from ' + postable + ' group by contig'
		table = cur.execute(sql).fetchall()
		chr2max=[]
		firstrow=True
		counttable=postable + '_' + reftable + '_count' + limit  
		counttable_cnt=postable + '_' + reftable + '_count' + limit 
		for row in table:
			contig=row[0]
			chrmax=row[1]
			print('counting ' + contig + ' ' + str(chrmax))
			rangestart=1
			rangeend=rangestart+byrange
			rangestartref=-int(limit)
			rangeendref=rangestart+byrange+int(limit)
			while rangestart<=chrmax:
				print(str(rangestart) + '-' + str(rangeend) + '  -> ' +   str(rangestartref) + '-' + str(rangeendref))
				sql='select pk.contig, pk.bp, count(rk.bp-pk.bp) as count_' + limit + ' from ' + postable + ' pk, ' + reftable + ' rk where pk.contig=rk.contig and abs(rk.bp-pk.bp)<=' + limit + " and rk.bp<>pk.bp and pk.contig='" + contig + "' and rk.contig='" + contig + "' and (pk.bp between " + str(rangestart) + ' and ' + str(rangeend) + ') and (rk.bp between ' + str(rangestartref) + ' and ' + str(rangeendref) + ') group by pk.contig,pk.bp;'
				if firstrow:
					sql='CREATE TABLE IF NOT EXISTS ' +counttable_cnt +  ' as ' + sql
					print(sql)
					if not printsqlonly:
						cur.execute(sql)				
						#cur.execute('CREATE TABLE IF NOT EXISTS ' + postable + '_' + reftable + '_count' + limit +  ' as select pk.contig, pk.bp, count(rk.bp-pk.bp) as count_' + limit + ' from ' + postable + ' pk, ' + reftable + ' rk where pk.contig=rk.contig and abs(rk.bp-pk.bp)<=' + limit + ' and rk.bp<>pk.bp group by pk.contig,pk.bp;')
						cur.execute("CREATE INDEX IF NOT EXISTS index_" + counttable_cnt  + "_contigbp ON " +  counttable_cnt   + "(contig, bp);")
						cur.execute("CREATE INDEX IF NOT EXISTS index_" + counttable_cnt  + "_contig ON " + counttable_cnt + "(contig);")
						cur.execute("CREATE INDEX IF NOT EXISTS index_" + counttable_cnt  + "_bp ON " + counttable_cnt + "(bp);")
					firstrow=False
				else:
					sql='INSERT INTO ' + counttable_cnt  + ' ' + sql
					print(sql)
					if not printsqlonly:
						cur.execute(sql)				
				rangestart=rangestart+byrange
				rangeend=rangeend+byrange
				rangestartref=rangestartref+byrange
				rangeendref=rangeendref+byrange

	else:
		create_counts_table(postable,reftable,limit)

	if withintable and counttable and withinlimit:
		withinlimit=str(withinlimit)
		withincounttable=withintable + '_' + counttable + '_within' + withinlimit
		sql='CREATE TABLE IF NOT EXISTS ' + withincounttable + ' as select pk.contig, pk.bp, count(rk.bp-pk.bp) as count_' + withinlimit + ' from ' + withintable + ' pk, ' + counttable + ' rk where pk.contig=rk.contig and abs(rk.bp-pk.bp)<=' + withinlimit + ' and rk.bp<>pk.bp group by pk.contig,pk.bp;'
		print(sql)
		if not printsqlonly:
			cur.execute(sql)
			cur.execute("CREATE INDEX IF NOT EXISTS index_" + postable + '_' + reftable + '_count' + limit + '_' + withintable  + "_contigbp ON " + postable + '_' + reftable + '_count' + limit  + "(contig, bp);")
			cur.execute("CREATE INDEX IF NOT EXISTS index_" + postable + '_' + reftable + '_count' + limit + '_' + withintable  + "_contig ON " + postable + '_' + reftable + '_count' + limit  + "(contig);")
			cur.execute("CREATE INDEX IF NOT EXISTS index_" + postable + '_' + reftable + '_count' + limit + '_' + withintable  + "_bp ON " + postable + '_' + reftable + '_count' + limit  + "(bp);")




	print('loaded ' + postable + '_' + reftable + '_count' + limit)	



def create_counts_summary_adj10(postable, reftable, withinlimit, adjlimit):
	# create the final result table with withinlimit, adjlimit
	counttable=reftable + '_' + reftable + '_count' + str(adjlimit) 
	posref=postable + '_' + reftable
	withintable=postable
	withincounttable=withintable + '_' + counttable + '_within' + str(withinlimit)

	sql='create table if not exists ' +  posref + '_flanks as SELECT pk.contig, pk.bp, (CASE WHEN pr10.count_' + str(withinlimit) + ' IS NOT NULL then pr10.count_' + str(withinlimit) + ' ELSE 0 END) AS within100_adj10 , (CASE WHEN pr100.count_' + str(withinlimit) + ' IS NOT NULL then pr100.count_' + str(withinlimit) + ' ELSE 0 END) AS within100 FROM ' + postable + ' pk LEFT JOIN ' + withincounttable + ' pr10 ON pr10.contig=pk.contig AND pr10.bp=pk.bp  LEFT JOIN ' + posref + '_count100 pr100 ON pr100.contig=pk.contig AND pr100.bp=pk.bp order by pk.contig, pk.bp'
	
	print(sql)
	if not printsqlonly:
		cur.execute(sql)
		cur.execute("CREATE INDEX IF NOT EXISTS index_" + posref + '_flanks_contigbp ON ' +  posref + '_flanks (contig,bp)') 
	print('loaded ' + posref + '_flanks')		
	return  posref + '_flanks'


def create_counts_summary(postable, reftable):
	# create the final result table if using only create_counts_table 
	posref=postable + '_' + reftable
	sql='create table if not exists ' + posref + '_flanks as SELECT pk.contig, pk.bp, (CASE WHEN pr10.count_10 IS NOT NULL then pr10.count_10 ELSE 0 END) AS within100_adj10 , (CASE WHEN pr100.count_100 IS NOT NULL then pr100.count_100 ELSE 0 END) AS within100 FROM ' + postable + ' pk LEFT JOIN ' + posref + '_count10 pr10 ON pr10.contig=pk.contig AND pr10.bp=pk.bp  LEFT JOIN ' + posref + '_count100 pr100 ON pr100.contig=pk.contig AND pr100.bp=pk.bp order by pk.contig, pk.bp'
	print(sql)
	cur.execute(sql)
	cur.execute("CREATE INDEX IF NOT EXISTS index_" + posref + '_flanks_contigbp ON ' +  posref + '_flanks (contig,bp)') 
	print('loaded ' + posref + '_flanks')		



create_table(posname,poslist,inrefcols=poscols)
create_table(refname,reflist,inrefcols=refcols)
create_counts_table_byrange(posname,refname,100,byrange)
create_counts_table_byrange(refname,refname,10,byrange,posname,100)
posrefadj_table= create_counts_summary_adj10(posname, refname, 100, 10)
con.commit()


if not printsqlonly:
	db_df = pd.read_sql_query("SELECT * FROM " + posrefadj_table + " order by contig,bp", con)

	outfle='countflankingsnps-100-10-' + posrefadj_table +  '.txt'
	if includechr:
		outfle='countflankingsnps-100-10-' + posrefadj_table + '-' + includechr + '.txt'
	if excludechr:
		outfle='countflankingsnps-100-10-' + posrefadj_table + '-' + excludechr + '.txt'

	db_df.to_csv(outfle, index=False, sep='\t')
	print('generated ' + outfle)

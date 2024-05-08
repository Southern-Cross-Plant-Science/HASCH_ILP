import subprocess
import os

# for large SNP list, this scripit will start the main program countflankingsnps.py by chromosome + all unassembled contigs 

# SNP list to count flanking regions
poslist='maf2_fmis6.isec_ipkgbs_5.vcf.gz.plinkA.map'
# all assembled chromosome name starts with
assembled_startswith="NC_"
# then edit the settings in countflankingsnps.py



output = subprocess.check_output("cut -f1 " + poslist + " | sort | uniq", shell=True).decode('utf-8')
output=output.split('\n')
output.remove("")

outputs=[]
for ctg in output:
	if ctg.startswith(assembled_startswith):
		cmd='nohup python -u count-flankings-11_clean_2.py ' + ctg + " > nohup-count-flankings-" + ctg + '.out &'
		print(cmd)
		os.system(cmd)
		outputs.append('countflankingsnps-100-10-pos57kc_ref5m_flanks-' + cgt + '.txt')

cmd='nohup python -u countflankingsnps.py "" ' + assembled_startswith + ' >  nohup-count-flankings-' + assembled_startswith + '.out &'
print(cmd)
os.system(cmd)		
outputs.append('countflankingsnps-100-10-pos57kc_ref5m_flanks-' + assembled_startswith + '.txt')

print('when all are finished, concatenate all outputs' )
print('cat ' + " ".join(outputs) + " > countflankingsnps-100-10-pos57k_ref5m_flanks-all.txt")


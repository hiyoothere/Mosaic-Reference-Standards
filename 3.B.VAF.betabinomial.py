import glob
import scipy
from scipy import stats
from scipy.stats import betabinom
import sys

#DataPath = sys.argv[1]
DataPath="/data/project/MRS/0.Genotype/4.analysis/2.NC/germline/4.filter_binomial_af/interim"
#OutPath = sys.argv[2]
OutPath="/data/project/MRS/0.Genotype/4.analysis/2.NC/germline/4.filter_binomial_af/interim"
filelist =  glob.glob(DataPath + '/germline.*het*.vcf')
#TAG = "het_ind"
print(filelist)

valid_out = open("NEW.valid.betabinom_74_76.het_ind.txt", 'w')
invalid_out = open("NEW.invalid.betabinom_74_76.het_ind.txt", 'w')
valid_out.write('\t'.join(["chr", "pos", "depth", "AF", "Beta_prob","p_value"])+'\n')
invalid_out.write('\t'.join(["chr", "pos", "depth", "AF", "Beta_prob", "p_value"])+'\n')


temp_ct = 0
for file in filelist:
	print(file)
	outname = "bB_0.01."+str(file).split('/')[-1]
    
	output = open(OutPath +'/' + outname , 'wb')	
	f = open(file,'rb')

	## count holders
	total = 0
	valid = 0
	invalid = 0

	for og_line in f:
		
		line = og_line.decode("utf-8")
		if ('CHROM' in line):
			pass
			#output.write(line)
		elif ("fileformat" in line):
			pass
			#output.write(line)
		elif line[0] != '#':
			info = line.split('\t')[-1]
			s = line.split()
			
			dp = int(s[2])
			af = float(s[5]) #jkim1105 : no alt ct, so using af*dp
			alt = int(dp * af) #jkim1105 : no alt ct, so using af*dp
			
			### cumulative Betabinomial distribution
			### alpha / beta of S0 is 74 and 76
			prob = betabinom.cdf(alt, dp, 74,76) 
			
			if prob < 0.5:
				p_value = prob*2
			else:
				p_value = (1 - prob)*2

			if p_value >= 0.01:
				if dp >= 40:
					valid += 1
					output.write(og_line)
					#valid_out.write('\t'.join([s[0], s[1], s[2], s[5], str(prob), str(p_value)]) + '\n')
			else:
				invalid += 1
				temp_ct +=1
				#invalid_out.write('\t'.join([s[0], s[1], s[2], s[5], str(prob), str(p_value)]) + '\n')
		else:
			print(line)
			pass
			#invalid += 1
			#invalid_out.write('\t'.join([s[0], s[1], s[2], s[5], str(p_value)]) + '\n')
		
		total += 1

	valid_out.close()
	invalid_out.close()
	
	output.close()
	f.close()
	print(temp_ct)

	print(str(file))
	print("=======RESULT=======")
	print("total is : " + str(total))
	print("valid is : " + str(valid))
	print("invalid is : " + str(invalid))
	print(valid + invalid)
 

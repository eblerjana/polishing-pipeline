import sys
import gzip

def extract_genotype_alleles(line):
	"""
	Extract genotype alleles.
	"""
	fields = line.strip().split()
	gt_index = fields[8].split(':').index("GT")
	genotype_str = fields[9].split(':')[gt_index]
	if '.' in genotype_str:
		return None,None
	alleles = [fields[3]] + fields[4].split(',')
	gt_alleles = [int(a) for a in genotype_str.replace('/', '|').split('|')]
	assert len(gt_alleles) == 2
	return (fields[0], fields[1]), sorted([alleles[gt_alleles[0]], alleles[gt_alleles[1]]])


vcf = sys.argv[1]
var_to_gt = {}

# read first VCF and store genotype alleles
for line in gzip.open(vcf, 'rt'):
	if line.startswith('#'):
		continue
	pos, gt_alleles = extract_genotype_alleles(line)	
	var_to_gt[pos] = gt_alleles


# read second VCF from stdin and keep only variants with matching genotypes
for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	pos, gt_alleles = extract_genotype_alleles(line)
	if pos is None:
		continue
	if not pos in var_to_gt:
		continue
	if var_to_gt[pos] != gt_alleles:
		sys.stderr.write("Skipping variant " + pos[0] + ":" + pos[1] + " with mismatching genotypes.\n")
		continue
	print(line.strip())

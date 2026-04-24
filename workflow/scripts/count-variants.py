import sys

nr_undefined = 0
nr_het = 0
nr_hom = 0
nr_absent = 0
nr_total = 0

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	gt_index = fields[8].split(':').index('GT')
	gt = fields[9].strip().split(':')[gt_index]

	nr_total += 1

	if '.' in gt:
		nr_undefined += 1
		continue

	alleles = gt.replace('|', '/').split('/')

	if (alleles[0] == "0") and (alleles[1] == "0"):
		nr_absent += 1
		continue

	if alleles[0] == alleles[1]:
		nr_hom += 1
		continue
	nr_het += 1

print('\t'.join(['nr_absent', 'nr_het', 'nr_hom', 'nr_undefined', 'nr_total']))
print('\t'.join([str(nr_absent), str(nr_het), str(nr_hom), str(nr_undefined),  str(nr_total)]))

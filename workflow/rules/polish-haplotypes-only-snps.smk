
rule polish_phase_variants_only_snps:
	"""
	Phase variants using WhatsHap and ONT reads.
	"""
	input: 
		chrom = "{results}/polishing/{callset}/contig-names/{sample}_{haplotype}/{chrom}.txt",
		bam = "{results}/polishing/{callset}/ont/ont_{sample}_{haplotype}_{chrom}.bam",
		vcf = "{results}/polishing/{callset}/deepvariant/deepvariant_{sample}_{haplotype}_{chrom}_filtered.vcf.gz",
		reference = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		vcf_gz = "{results}/polishing/{callset}/only-snps/whatshap/whatshap_{sample}_{haplotype}_{chrom}.vcf.gz",
		vcf = temp("{results}/polishing/{callset}/only-snps/whatshap/whatshap_{sample}_{haplotype}_{chrom}.vcf")
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	resources:
		mem_mb = 50000,
		walltime = "04:00:00"
	threads: 1
	log:
		"{results}/polishing/{callset}/only-snps/whatshap/whatshap_{sample}_{haplotype}_{chrom}.log"
	benchmark:
		"{results}/polishing/{callset}/only-snps/whatshap/whatshap_{sample}_{haplotype}_{chrom}.benchmark.txt"
	shell:
		" whatshap phase --only-snvs --chromosome {wildcards.chrom} -o {output.vcf} --reference={input.reference} {input.vcf} {input.bam} &> {log} "
		" && "
		" bgzip -c {output.vcf} > {output.vcf_gz} "
		" && "
		" tabix -p vcf {output.vcf_gz} "


rule polish_synchronize_variants_only_snps:
	"""
	Synchronize phased blocks and filter VCF for haplotagging.
	"""
	input:
		"{results}/polishing/{callset}/only-snps/whatshap/whatshap_{sample}_{haplotype}_{chrom}.vcf.gz"
	output:
		"{results}/polishing/{callset}/only-snps/haplotag/haplotag_{sample}_{haplotype}_{chrom}.vcf.gz"
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	benchmark:
		"{results}/polishing/{callset}/only-snps/haplotag/haplotag_{sample}_{haplotype}_{chrom}.benchmark.txt"
	log:
		"{results}/polishing/{callset}/only-snps/haplotag/haplotag_{sample}_{haplotype}_{chrom}.log"
	shell:
		"""
		python3 workflow/scripts/synchronize-phased-blocks.py {input} 2> {log} | bgzip > {output}
		tabix -p vcf {output}
		"""


rule polish_detect_errors_only_snps:
	"""
	Detect errors in the consensus haplotypes.
	"""
	input:
		small = "{results}/polishing/{callset}/only-snps/haplotag/haplotag_{sample}_{haplotype}_{chrom}.vcf.gz"
	output:
		vcf = "{results}/polishing/{callset}/only-snps/errors/errors_{sample}_{haplotype}_{chrom}_{windowsize}.vcf.gz"
	wildcard_constraints:
		haplotype = "hap1|hap2"
	conda:
		"../envs/whatshap.yml"
	log:
		"{results}/polishing/{callset}/only-snps/errors/errors_{sample}_{haplotype}_{chrom}_{windowsize}.log"
	benchmark:
		"{results}/polishing/{callset}/only-snps/errors/errors_{sample}_{haplotype}_{chrom}_{windowsize}.benchmark.txt"
	shell:
		" bcftools view -v snps {input.small} | python3 workflow/scripts/find_errors.py -window-width {wildcards.windowsize} 2> {log} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output.vcf} "
		" && "
		" tabix -p vcf {output.vcf} "




def aggregate_concat_phased_variants_only_snps(wildcards):
	checkpoint_output = checkpoints.polish_determine_contig_names.get(**wildcards).output[0]
	return expand("{results}/polishing/{callset}/only-snps/errors/errors_{sample}_{haplotype}_{chrom}_{windowsize}.vcf.gz",
			results=wildcards.results,
			callset=wildcards.callset,
			sample=wildcards.sample,
			haplotype=wildcards.haplotype,
			chrom=glob_wildcards(os.path.join(checkpoint_output, "{chrom}.txt")).chrom,
			windowsize=wildcards.windowsize)


rule polish_concat_phased_variants_only_snps:
	"""
	Concat chromosome-wise phased VCFs into a single VCF.
	"""
	input:
		vcfs = aggregate_concat_phased_variants_only_snps
	output:
		"{results}/polishing/{callset}/only-snps/errors/errors_{sample}_{haplotype}_{windowsize}.vcf.gz"
	conda:
		"../envs/whatshap.yml"
	wildcard_constraints:
		haplotype = "hap1|hap2",
		windowsize = "|".join([str(w) for w in WINDOW_WIDTHS])
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	threads: 15
	log:
		"{results}/polishing/{callset}/only-snps/errors/errors_whatshap_{sample}_{haplotype}_{windowsize}.log"
	benchmark:
		"{results}/polishing/{callset}/only-snps/errors/errors_whatshap_{sample}_{haplotype}_{windowsize}.benchmark.txt"
	shell:
		"""
		bcftools concat -o {output} -Oz --threads {threads} {input.vcfs} &> {log}
		bcftools index {output}
		"""

rule polish_correct_errors_only_snps:
	"""
	Correct errors of the consensus haplotypes.
	"""
	input:
		errors = "{results}/polishing/{callset}/only-snps/errors/errors_{sample}_{haplotype}_{windowsize}.vcf.gz",
		haplotype = "{results}/polishing/{callset}/consensus/{sample}_{haplotype}.fa"
	output:
		temp("{results}/polishing/{callset}/only-snps/polished/{sample}_{haplotype}_{windowsize}.fa")
	log:
		"{results}/polishing/{callset}/only-snps/polished/{sample}_{haplotype}_{windowsize}.log"
	benchmark:
		"{results}/polishing/{callset}/only-snps/polished/{sample}_{haplotype}_{windowsize}.benchmark.txt"
	conda:
		"../envs/whatshap.yml"
	shell:
		"""
		bcftools consensus --haplotype A -f {input.haplotype}  -e \'ALT~\"<.*>\"\' {input.errors} 2> {log}  > {output}
		"""

rule polish_compress_haplotypes_only_snps:
	input:
		expand("{{results}}/polishing/{{callset}}/only-snps/polished/{sample}_{haplotype}_{{windowsize}}.fa", sample = SAMPLES, haplotype = ["hap1", "hap2"])
	output:
		"{results}/polishing/{callset}/only-snps/polished/{callset}_{windowsize}_polished.agc"
	benchmark:
		"{results}/polishing/{callset}/only-snps/polished/{callset}_{windowsize}_polished.benchmark.txt"
	log:
		"{results}/polishing/{callset}/only-snps/polished/{callset}_{windowsize}_polished.log"
	conda:
		"../envs/bwa.yml"
	threads:
		24
	resources:
		mem_mb = 100000,
		walltime = "02:00:00"
	shell:
		"""
		agc create {input} -o {output} -t {threads} &> {log}
		"""





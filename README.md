# Polishing Pipeline

Pipeline to polish consensus haplotypes generated from phased genotypes based on short and long read data.

# Input data

* **consensus haplotypes**: AGC archive containing FASTA sequences of haplotypes to be polished.
* **haplotype names**: TSV file specifying sequence names of the haplotypes in the AGC archive. If not provided, sequence names are assumed to follow the pattern: `` <callset>_<sample>_<H1|H2> ``. Required columns (in that order):

  ```
   <callset name>     <sample>_<H1|H2>      <sequence name in AGC>
  ```

* **sample sheet**: TSV file providing paths to Illumina and long read (ONT for HiFi) data. Required columns (in that order):
  ```
   <sample name> <path/to/illumina.fasta/fastq> <path/to/long-read.cram> <read tech: ONT|HIFI>
  ```
  Currently, short read are expected in FASTA/FASTQ format, while long reads (ONT or HiFi) are expected in CRAM format.
* **cram reference**: Reference FASTA file with reference sequence used for long read cram.
* **list of haploid chromosomes**: File providing a list of all chromosomes in the consensus haplotypes that shall be treated as haploid. Required columns:

  ```
  <sample>        <comma sep. names of haploid chroms>
  ```

# How to run

Provide paths to the input files in ``config/config.yaml``:

```
# name of the output folder
results: "results"

# archive with compressed haplotypes to be polished
agc:
   CHRY: "path/to/consensus.agc"

# tab-delimited file: <trio name>	<child>	<father>	<mother>	<sex>	<population>	<superpopulation>	<path to Illumina FASTA>	<path to long read CRAM>	<long read tech ONT|HIFI>
sample_sheet: "path/to/sample-sheet.tsv"

# reference used in long read CRAMs
cram_ref: "path/to/cram_ref.fa"

# optional: tab-delimited file: <sample>        <comma sep. names of haploid chroms>
# if not available, set to ""
haploid_chroms: "path/to/haploid-chroms.txt"

# tab-delimited file of format: <callset>     <sample>_<H1|H2>      name_in_agc
# specifying the name of each haplotype sequence in the AGC archive.
# if not provided, pipeline assumes names follow this pattern: <callset>_<sample>_<H1|H2>
haplotype_names: "path/to/names.tsv"
```

Then, run the pipeline using snakemake, e.g. using:

```
snakemake --use-conda --use-singularity -j <nr_cores>
```

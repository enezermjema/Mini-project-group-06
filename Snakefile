SAMPLES, = glob_wildcards("Data/{sample}.fastq")

Illumina_ref = "Data/Lablab_ref/Lablab_purpureus_147D_AOCC.fa"
Nanopore_ref = "Data/Lablab_ref/isaac_147D_flye2.7_hypo-polished_assembly_2020.fa"

# wildcard_constraints:
#    seqfile = '\w+'

rule all:
	input:
		expand("Analyses/variants/illumina/{sample}_stats.txt", sample=SAMPLES),
# 		expand("Analyses/Illumina_variants/SV_sniff/{sample}.vcf_stats.txt", sample=SAMPLES)

rule nanoqc:      # for the 147X 
	input: 
		"Data/{sample}.fastq"
	output: 
		"Analyses/nanoqc_report/{sample}",
		touch("Analyses/nanoqc_logs/{sample}_NanoQC.log")
	shell: 
		"nanoQC -o {output[0]} {input}"


rule porechop:           # trimming adapters for 147X
	input:
		"Analyses/nanoqc_logs/{sample}_NanoQC.log",
		"Data/{sample}.fastq"
	output:
		"Analyses/trimmed_reads/{sample}_trimmed.fastq"
	threads: 4
	shell:
		"porechop -i {input[1]} -t {threads} -o {output}"

rule minimap2_illum_nanop:       # index 147D illumina ref and map 147X reads to illumina 147Dref
	input:
		ref_illum = Illumina_ref,
		ref_nanop = Nanopore_ref,
		READS = "Analyses/trimmed_reads/{sample}_trimmed.fastq"
	output:
		"Analyses/minimap2_mapping/illumina/{sample}_mapped.sam",
		"Analyses/minimap2_mapping/nanopore/{sample}_mapped.sam"
	threads: 4
	shell:
		"""
		minimap2 -t {threads} -ax map-ont {input.ref_illum} {input.READS} > {output[0]}
		minimap2 -t {threads} -ax map-ont {input.ref_nanop} {input.READS} > {output[1]}
		"""

	
rule sam_sorting1:              # samtools sort and samtools index  -> sam file to indexed bam
	input:
		ref_illum = Illumina_ref,
		ref_nanop = Nanopore_ref,
		illum_file = "Analyses/minimap2_mapping/illumina/{sample}_mapped.sam",
		nanop_file = "Analyses/minimap2_mapping/nanopore/{sample}_mapped.sam"
	output:
		"Analyses/minimap2_mapping/illumina/{sample}_mapped_sorted.bam",
		"Analyses/minimap2_mapping/nanopore/{sample}_mapped_sorted.bam"
	threads: 4
	shell:
		"""
		samtools sort -t {threads} -o {output[0]} --reference {input.ref_illum} {input.illum_file}
		samtools index {output[0]}

		samtools sort -t {threads} -o {output[1]} --reference {input.ref_nanop} {input.nanop_file}
		samtools index {output[1]}
		"""

# variant calling (SNPs and indels)
rule freebayes_vc:       
	input:
		ref_illum = Illumina_ref,
		ref_nanop = Nanopore_ref,
		illum_bam = "Analyses/minimap2_mapping/illumina/{sample}_mapped_sorted.bam",
		nanop_bam = "Analyses/minimap2_mapping/nanopore/{sample}_mapped_sorted.bam"
	output:
		"Analyses/variants/illumina/{sample}.vcf",
		"Analyses/variants/nanopore/{sample}.vcf"
	shell:
		"""
		freebayes -f {input.ref_illum} {input.illum_bam} > {output[0]}
		freebayes -f {input.ref_nanop} {input.nanop_bam} > {output[1]}
		
		"""

rule vc_metrics:
	input:
		illum_vcf = "Analyses/variants/illumina/{sample}.vcf",
		nanop_vcf = "Analyses/variants/nanopore/{sample}.vcf"
	output:
		"Analyses/variants/illumina/{sample}_stats.txt",
		"Analyses/variants/nanopore/{sample}_stats.txt"
	shell:
		"""
		bgzip {input.illum_vcf} && tabix -p vcf {input.illum_vcf}.gz | rtg vcfstats > {output[0]}
		bgzip {input.nanop_vcf} && tabix -p vcf {input.nanop_vcf}.gz | rtg vcfstats > {output[1]}

		"""

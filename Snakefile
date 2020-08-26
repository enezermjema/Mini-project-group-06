SAMPLES, = glob_wildcards("Data/{sample}.fastq")

Illumina_ref = "Data/Lablab_ref/Lablab_purpureus_147D_AOCC.fa"
Nanopore_ref = "Data/Lablab_ref/isaac_147D_flye2.7_hypo-polished_assembly_2020.fa"

# wildcard_constraints:
#    seqfile = '\

rule all:
	input:
		expand("Analyses/variants/illumina/{sample}_stats.txt", sample=SAMPLES),
		expand("Analyses/variants/illumina/structural_variants/{sample}_stats.txt", sample=SAMPLES)

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




# structural variant calling

# mapping with ngmlr for Svc
rule ngmlr_svc:
	input: 
		ref_illum = Illumina_ref,
		ref_nanop = Nanopore_ref,
		READS = "Analyses/trimmed_reads/{sample}_trimmed.fastq"
	output: 
		"Analyses/ngmlr_mapping/illumina/{sample}_ngmlr_mapped.sam",
		"Analyses/ngmlr_mapping/nanopore/{sample}_ngmlr_mapped.sam"
	threads: 4
	shell: 
		"""
		ngmlr -t {threads} -r {input.ref_illum}  -q {input.READS} -o {output[0]} -x ont
		ngmlr -t {threads} -r {input.ref_nanop}  -q {input.READS} -o {output[1]} -x ont
		"""


rule sam_sorting2:              # samtools sort and samtools index  -> sam file to indexed bam
	input:
		ref_illum = Illumina_ref,
		ref_nanop = Nanopore_ref,
		illum_fil = "Analyses/ngmlr_mapping/illumina/{sample}_ngmlr_mapped.sam",
		nanop_fil = "Analyses/ngmlr_mapping/nanopore/{sample}_ngmlr_mapped.sam"
	output:
		"Analyses/ngmlr_mapping/illumina/{sample}_ngmlr_mapped_sorted.bam",
		"Analyses/ngmlr_mapping/nanopore/{sample}_ngmlr_mapped_sorted.bam"
	threads: 4
	shell:
		"""
		samtools sort -t {threads} -o {output[0]} --reference {input.ref_illum} {input.illum_fil}
		samtools index {output[0]}

		samtools sort -t {threads} -o {output[1]} --reference {input.ref_nanop} {input.nanop_fil}
		samtools index {output[1]}
		"""

rule sniffles_svc:
	input: 
		illum_ba = "Analyses/ngmlr_mapping/illumina/{sample}_ngmlr_mapped_sorted.bam",
		nanop_ba = "Analyses/ngmlr_mapping/nanopore/{sample}_ngmlr_mapped_sorted.bam"
	output: 
		"Analyses/variants/illumina/structural_variants/{sample}_sv.vcf",
		"Analyses/variants/nanopore/structural_variants/{sample}_sv.vcf"
	shell: 
		"""
		sniffles  -m {input.illum_ba} -v {output[0]}
		sniffles  -m {input.nanop_ba} -v {output[1]}
		"""


rule sniff_sv_met:
	input:
		illum_vcf1 = "Analyses/variants/illumina/structural_variants/{sample}_sv.vcf",
		nanop_vcf1 = "Analyses/variants/nanopore/structural_variants/{sample}_sv.vcf"
	output:
		"Analyses/variants/illumina/structural_variants/{sample}_stats.txt",
		"Analyses/variants/nanopore/structural_variants/{sample}_stats.txt"
	shell:
		"""

		bgzip {input.illum_vcf1} && tabix -p vcf {input.illum_vcf1}.gz | rtg vcfstats > {output[0]}
		bgzip {input.nanop_vcf1} && tabix -p vcf {input.nanop_vcf1}.gz | rtg vcfstats > {output[1]}	

		"""



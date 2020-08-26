
SAMPLES, = glob_wildcards("Data/{sample}.fastq")

Illumina_ref = "Data/Lablab_ref/Lablab_purpureus_147D_AOCC.fa"
Nanopore_ref = "Data/Lablab_ref/isaac_147D_flye2.7_hypo-polished_assembly_2020.fa"

<<<<<<< HEAD
=======
# wildcard_constraints:
#    seqfile = '\
>>>>>>> a897bce9548d2a273e546f9ee9c04e47fe415e5f

rule all:
	input:
		expand("Analyses/variants/illumina/{sample}_stats.txt", sample=SAMPLES),
		expand("Analyses/variants/illumina/structural_variants/{sample}_stats.txt", sample=SAMPLES)

rule nanoqc:      # for the 147X 
	input: 
		"Data/{sample}.fastq"
	output: 
		directory("Analyses/nanoqc_report/{sample}"),
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
<<<<<<< HEAD
=======




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
>>>>>>> a897bce9548d2a273e546f9ee9c04e47fe415e5f
		sniffles  -m {input.illum_ba} -v {output[0]}
		sniffles  -m {input.nanop_ba} -v {output[1]}
		"""


<<<<<<< HEAD
# generating variants statistics 
rule variants_stats:
    input:
        illum_vcf = "Analyses/variants/illumina/{sample}.vcf",
        nanop_vcf = "Analyses/variants/nanopore/{sample}.vcf",
        illum_sv = "Analyses/variants/illumina/structural_variants/{sample}_sv.vcf",
        nanop_sv = "Analyses/variants/nanopore/structural_variants/{sample}_sv.vcf"
    output:
        "Analyses/variants/illumina/{sample}_stats.txt",
        "Analyses/variants/nanopore/{sample}_stats.txt",
        "Analyses/variants/illumina/structural_variants/{sample}_stats.txt",
        "Analyses/variants/nanopore/structural_variants/{sample}_stats.txt"
    shell:
        """
        rtg vcfstats {input[0]} > {output[0]}
        rtg vcfstats {input[1]} > {output[1]}
        rtg vcfstats {input[2]} > {output[2]}
        rtg vcfstats {input[3]} > {output[3]}  
        """


=======
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


=======
SAMPLES, = glob_wildcards("Data/{samples}.fastq")

rule all:
    input:
        expand("Analyses/Variant_calling/sniffles_illumina_{samples}.vcf", samples=SAMPLES),
        expand("Analyses/Variant_calling/sniffles_nanopore_{samples}.vcf", samples=SAMPLES)

rule nanoqc:
    input:
        "Data/{samples}.fastq"
    output: 
        "{Analyses}/{nanopore_QC}/{samples}.fastq.html",                
    shell:
        "nanoQC -o {output} {input}"

rule porechop:
    input:
        "Data/{samples}.fastq"
    output:
        "Analyses/nanopore_QC/trimmed_{samples}.fastq"
    threads: 6
    shell:
        "porechop -i {input} -t {threads} -o {output}"

rule reference_indexing:
    input:
        Illumina_ref = "Data/Lablab_ref/Lablab_purpureus_147D_AOCC.fa",
        Nanopore_ref = "Data/Lablab_ref/isaac_147D_flye2.7_hypo-polished_assembly_2020.fa"
    output:
        Illumina_fai = "Data/Lablab_ref/Lablab_purpureus_147D_AOCC.fa.fai",
        Nanopore_fai = "Data/Lablab_ref/isaac_147D_flye2.7_hypo-polished_assembly_2020.fa.fai"
    shell:
        "samtools faidx {input} > {output}"
#       "samtools faidx {input.Illumina_ref} > {output.Illumina_fai}"
#        "samtools faidx {input.Nanopore_ref} > {output.Nanopore_fai}"

rule alignment:
    input:
        Illumina_ref = "Data/Lablab_ref/Lablab_purpureus_147D_AOCC.fa",
        Nanopore_ref = "Data/Lablab_ref/isaac_147D_flye2.7_hypo-polished_assembly_2020.fa",
        reads = "Analyses/nanopore_QC/trimmed_{samples}.fastq"
    output:
        Illumina_sam = "Analyses/{Alignment}/mapped_illumina_{samples}.sam",
        Nanopore_sam = "Analyses/{Alignment}/mapped_nanopore_{samples}.sam"
    threads: 8
    shell:
        "minimap2 -t {threads} -ax map-ont {input.Illumina_ref} {input.reads} > {output.Illumina_sam}; minimap2 -t {threads} -ax map-ont {input.Nanopore_ref} {input.reads} > {output.Nanopore_sam}"
#        "minimap2 -t {threads} -ax map-ont {input.Nanopore_ref} {input.reads} > {output.Nanopore_sam}"

rule sam_sorting:
    input:
        Illumina_sam = "Analyses/Alignment/mapped_illumina_{samples}.sam",
        Nanopore_sam = "Analyses/Alignment/mapped_nanopore_{samples}.sam"
    output:
        Illumina_bam = "Analyses/Alignment/sorted_illumina_{samples}.bam",
        Nanopore_bam = "Analyses/Alignment/sorted_nanopore_{samples}.bam"
    shell:
        "samtools view -O BAM {input.Illumina_sam} | samtools sort -O BAM -o {output.Illumina_bam}; samtools view -O BAM {input.Nanopore_sam} | samtools sort -O BAM -o {output.Nanopore_bam}"
#        "samtools view -O BAM {input.Nanopore_sam} | samtools sort -O BAM -o {output.Nanopore_bam}"

rule bam_indexing:
    input:
        Illumina_bam = "Analyses/Alignment/sorted_illumina_{samples}.bam",
        Nanopore_bam = "Analyses/Alignment/sorted_nanopore_{samples}.bam"
    output:
        Illumina_idx = "Analyses/Alignment/sorted_illumina_{samples}.bam.bai",
        Nanopore_idx = "Analyses/Alignment/sorted_Nanopore_{samples}.bam.bai"
    shell:
        "samtools index {input.Illumina_bam} > {output.Illumina_idx}; samtools index {input.Nanopore_bam} > {output.Nanopore_idx}"
#        "samtools index {input.Nanopore_bam}"

rule variant_calling:
    input:
        Illumina_bam = "Analyses/Alignment/sorted_illumina_{samples}.bam",
        Nanopore_bam = "Analyses/Alignment/sorted_nanopore_{samples}.bam",
        Illumina_ref = "Data/Lablab_ref/Lablab_purpureus_147D_AOCC.fa",
        Nanopore_ref = "Data/Lablab_ref/isaac_147D_flye2.7_hypo-polished_assembly_2020.fa"
    output:
        Illumina_vcf = "Analyses/{Variant_calling}/freebayes_illumina_{samples}.vcf",
        Nanopore_vcf = "Analyses/{Variant_calling}/freebayes_nanopore_{samples}.vcf"
    shell:
        "freebayes -f {input.Illumina_ref} {Illumina_bam} > {output.Illumina_vcf}; freebayes -f {input.Nanopore_ref} {Nanopore_bam} > {output.Nanopore_vcf}"
#        "freebayes -f {input.Nanopore_ref} {Nanopore_bam} > {output.Nanopore_vcf}"

rule StructuralVariant:
    input:
        Illumina_bam = "Analyses/Alignment/sorted_illumina_{samples}.bam",
        Nanopore_bam = "Analyses/Alignment/sorted_nanopore_{samples}.bam"
    output:
        Illumina_SV = "Analyses/Variant_calling/sniffles_illumina_{samples}.vcf",
        Nanopore_SV = "Analyses/Variant_calling/sniffles_nanopore_{samples}.vcf"
    shell:
        "sniffles -m {input.Illumina_bam} -v {output.Illumina_SV}; sniffles -m {input.Nanopore_bam} -v {output.Nanopore_SV}"
>>>>>>> a897bce9548d2a273e546f9ee9c04e47fe415e5f

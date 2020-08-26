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

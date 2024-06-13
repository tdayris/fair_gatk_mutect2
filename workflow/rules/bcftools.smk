rule fair_gatk_mutect2_bcftools_norm_split_multiallelic:
    input:
        "tmp/fair_gatk_mutect2_snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.vcf",
        ref=lambda wildcards: select_fasta(wildcards),
    output:
        temp(
            "tmp/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3_000,
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_bcftools_norm_split_multiallelic",
            default="--multiallelics -any --rm-dup none",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/norm"


rule fair_gatk_mutect2_bcftools_view:
    input:
        "tmp/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2i_bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_bcftools_view", default="--with-header"
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/view"


rule fair_gatk_mutect2_bcftools_index:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz",
    output:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_bcftools_index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/index"


rule fair_gatk_mutect2_bcftools_mutect2_stats:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz",
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz.tbi",
    output:
        temp(
            "tmp/fair_gatk_mutect2_bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.stats.txt"
        ),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 7) * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        lookup_config(dpath="params/fair_gatk_mutect2_bcftools_stats", default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/stats"

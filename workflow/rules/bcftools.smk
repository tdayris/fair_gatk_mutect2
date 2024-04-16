rule fair_gatk_mutect2_bcftools_view:
    input:
        "tmp/fair_gatk_mutect2/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.vcf",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2/bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2/bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(dpath="params/fair_gatk_mutect2/bcftools/view", default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/view"


rule fair_gatk_mutect2_bcftools_index:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz",
    output:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2/bcftools_index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2/bcftools_index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/index"


rule fair_gatk_mutect2_bcftools_mutect2_stats:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz",
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz.tbi",
    output:
        temp(
            "tmp/fair_gatk_mutect2/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.stats.txt"
        ),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 7) * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        lookup_config(dpath="params/fair_gatk_mutect2/bcftools/stats", default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/stats"

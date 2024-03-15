rule fair_gatk_mutect_germline_bcftools_view:
    input:
        "tmp/fair_gatk_mutect_germline/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.vcf",
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
        "logs/fair_gatk_mutect_germline/bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(
            dpath="params/fair_gatk_mutect_germline/bcftools/view", within=config
        ),
    wrapper:
        "v3.5.0/bio/bcftools/view"


rule fair_gatk_mutect_germline_bcftools_index:
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
        "logs/fair_gatk_mutect_germline/bcftools_index/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/bcftools_index/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    wrapper:
        "v3.5.0/bio/bcftools/index"


rule fair_gatk_mutect_germline_bcftools_mutect2_stats:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz",
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz.tbi",
    output:
        temp(
            "tmp/fair_gatk_mutect_germline/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.stats.txt"
        ),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 7) * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect_germline/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        lookup(dpath="params/fair_gatk_mutect_germline/bcftools/stats", within=config),
    wrapper:
        "v3.5.0/bio/bcftools/stats"

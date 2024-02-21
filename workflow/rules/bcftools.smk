rule fair_gatk_mutect_germline_bcftools_mutect2_stats:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/{sample}.vcf.gz",
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/{sample}.vcf.gz.tbi",
    output:
        temp(
            "tmp/fair_gatk_mutect_germline/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.stats.txt"
        ),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 7) * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_germline/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        lookup(dpath="params/bcftools/stats", within=config),
    wrapper:
        "v3.3.6/bio/bcftools/stats"

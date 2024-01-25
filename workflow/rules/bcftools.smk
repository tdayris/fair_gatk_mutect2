rule bcftools_mutect2_stats:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz",
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz.tbi",
    output:
        temp(
            "tmp/bcftools/stats/mutect2/germline/{species}.{build}.{release}.{datatype}/{sample}.stats.txt"
        ),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 7) * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir="tmp",
    log:
        "logs/bcftools/stats/{species}.{build}.{release}.{datatype}/{sample}.gatk.mutect2.germline.log",
    benchmark:
        "benchmark/bcftools/stats/{species}.{build}.{release}.{datatype}/{sample}.gatk.mutect2.germline.tsv"
    params:
        config.get("params", {}).get("bcftools", {}).get("stats", ""),
    wrapper:
        "v3.3.3/bio/bcftools/stats"

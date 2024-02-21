rule fair_gatk_mutect_germline_gatk_germline_varianteval:
    input:
        unpack(get_gatk_germline_varianteval_input),
    output:
        temp(
            "tmp/fair_gatk_mutect_germline/gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.grp"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/fair_gatk_mutect_germline/gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {}).get("gatk", {}).get("varianteval", ""),
    wrapper:
        "v3.3.6/bio/gatk/varianteval"

rule fair_gatk_mutect2_gatk_germline_varianteval:
    input:
        unpack(get_gatk_germline_varianteval_input),
    output:
        temp(
            "tmp/fair_gatk_mutect2_gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.grp"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_gatk_varianteval",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/gatk/varianteval"

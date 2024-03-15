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
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect_germline/gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(
            dpath="params/fair_gatk_mutect_germline/gatk/varianteval", within=config
        ),
    wrapper:
        "v3.5.0/bio/gatk/varianteval"

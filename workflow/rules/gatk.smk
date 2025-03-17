"""
## Memory
Requires a job with at most 37463.47  Mb,
 on average 28051.73 ± 17300.02 Mb, 
on Gustave Roussy's HPC Flamingo, on a 21.0  Mb dataset.
## Time
A job took 0:01:28 to proceed,
on average 0:00:30 ± 0:00:24
"""


rule fair_gatk_mutect2_gatk_germline_varianteval:
    input:
        unpack(get_gatk_germline_varianteval_input),
    output:
        temp(
            "tmp/fair_gatk_mutect2_gatk_germline_varianteval/{species}.{build}.{release}.{datatype}/{sample}.grp"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 15_000 + 30_000,
        runtime=lambda wildcards, attempt: attempt * 30,
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
        "v5.8.3/bio/gatk/varianteval"

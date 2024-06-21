"""
Reported on Flamingo on ~150 datasets
* time 27s ± 40s
*  mem 5.9Go ± 3.5Go
"""


rule fair_gatk_mutect2_snpsift_vartype:
    input:
        vcf="tmp/fair_gatk_mutect2_snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.vcf",
    output:
        vcf=temp(
            "tmp/fair_gatk_mutect2_snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.vcf"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4_000 + 6_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_snpsift_vartype", default=""
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/snpsift/varType"

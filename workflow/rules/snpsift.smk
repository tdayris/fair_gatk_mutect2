rule fair_gatk_mutect2_snpsift_vartype:
    input:
        vcf="tmp/fair_gatk_mutect2/snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.vcf",
    output:
        vcf=temp(
            "tmp/fair_gatk_mutect2/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.vcf"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2/snpsift/vartype", default=""
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/snpsift/varType"

rule fair_gatk_mutect_germline_snpsift_vartype:
    input:
        vcf="tmp/fair_gatk_mutect_germline/snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.vcf",
    output:
        vcf=temp(
            "tmp/fair_gatk_mutect_germline/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.vcf"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect_germline/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/snpsift_vartype/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(
            dpath="params/fair_gatk_mutect_germline/snpsift/vartype", within=config
        ),
    wrapper:
        "v3.5.0/bio/snpsift/varType"

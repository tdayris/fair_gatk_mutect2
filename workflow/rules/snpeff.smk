rule fair_gatk_mutect_germline_snpeff_download_reference:
    output:
        directory("reference/{species}.{build}.{release}/snpeff/{build}.{release}"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect_germline/snpeff_download_reference/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/snpeff_download_reference/{species}.{build}.{release}.tsv"
    params:
        reference="{build}.{release}",
    wrapper:
        "v3.5.0/bio/snpeff/download"


rule fair_gatk_mutect_germline_snpeff_annotate:
    input:
        calls="tmp/fair_gatk_mutect_germline_filter_mutect_calls/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz",
        calls_tbi="tmp/fair_gatk_mutect_germline_filter_mutect_calls/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz.tbi",
        db=getattr(
            lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=genomes,
            ),
            "snpeff_database",
            "reference/{species}.{build}.{release}/snpeff/{build}.{release}",
        ),
    output:
        calls=temp(
            "tmp/fair_gatk_mutect_germline/snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.vcf"
        ),
        stats="results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/QC/{sample}.html",
        csvstats=temp(
            "tmp/fair_gatk_mutect_germline/snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.csv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect_germline/snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup(dpath="params/fair_gatk_mutect_germline/snpeff", within=config),
    wrapper:
        "v3.5.0/bio/snpeff/annotate"

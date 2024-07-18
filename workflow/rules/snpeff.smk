"""
## Memory
Requires a job with at most 9927.27  Mb,
 on average 8509.23 ± 3751.78 Mb, 
on Gustave Roussy's HPC Flamingo, on a 1.0  Mb dataset.
## Time
A job took 0:18:33 to proceed,
on average 0:15:54 ± 0:07:00
"""


rule fair_gatk_mutect2_snpeff_download_reference:
    output:
        directory("reference/{species}.{build}.{release}/snpeff/{build}.{release}"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10_000,
        runtime=lambda wildcards, attempt: attempt * 60 * 3,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_snpeff_download_reference/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_snpeff_download_reference/{species}.{build}.{release}.tsv"
    params:
        reference="{build}.{release}",
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/snpeff/download"


"""
## Memory
Requires a job with at most 12370.43  Mb,
 on average 9273.7 ± 5713.04 Mb, 
on Gustave Roussy's HPC Flamingo, on a 21.0  Mb dataset.
## Time
A job took 0:09:39 to proceed,
on average 0:04:36 ± 0:03:02
"""


rule fair_gatk_mutect2_snpeff_annotate:
    input:
        calls="tmp/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz",
        calls_tbi="tmp/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz.tbi",
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
            "tmp/fair_gatk_mutect2_snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.vcf"
        ),
        stats="results/{species}.{build}.{release}.{datatype}/VariantCalling/QC/{sample}.html",
        csvstats=temp(
            "tmp/fair_gatk_mutect2_snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.csv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3_000 + 10_000,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(dpath="params/fair_gatk_mutect2_snpeff", default="-nodownload"),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/snpeff/annotate"

"""
Reported on Flamingo (hg38)
* time 11m (depends on bandwidth)
* mem 3.7Go
"""


rule fair_gatk_mutect2_snpeff_download_reference:
    output:
        directory("reference/{species}.{build}.{release}/snpeff/{build}.{release}"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4_000,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
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
Reported on Flamingo on 150 samples
* time 1:20   ± 10min
*  mem 8101mb ± 2Go
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
        mem_mb=lambda wildcards, attempt: attempt * 3_000 + 5_000,
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

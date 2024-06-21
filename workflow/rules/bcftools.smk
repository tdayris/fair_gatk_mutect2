"""
Reported on Flamingo on 150 datasets
* time 37s ± 31s
* mem 423mb ± 27mb
"""


rule fair_gatk_mutect2_bcftools_norm_split_multiallelic:
    input:
        "tmp/fair_gatk_mutect2_filter_mutect_calls/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz",
        ref=lambda wildcards: select_fasta(wildcards),
    output:
        temp(
            "tmp/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz",
        ),
        temp(
            "tmp/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz.tbi",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 200) + 500,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_norm_split_multiallelic/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_bcftools_norm_split_multiallelic",
            default="--multiallelics -any --rm-dup none --write-index='tbi'",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/norm"


"""
Reported on Flamingo on 150 datasets
* time  
* mem 
"""


rule fair_gatk_mutect2_bcftools_filter_pass:
    input:
        "tmp/fair_gatk_mutect2_snpeff_annotate/{species}.{build}.{release}.{datatype}/{sample}.vcf",
    output:
        temp(
            "tmp/fair_gatk_mutect2_bcftools_filter_pass/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 200 + 500,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_bcftools_filter_pass/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_filter_pass/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra='--include \'FILTER="PASS" || FILTER="."\'',
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/filter"


"""
Reported on Flamingo on 150 datasets
* time 13s ± 10s
* time 415mb ± 38mb
"""


rule fair_gatk_mutect2_bcftools_view:
    input:
        "tmp/fair_gatk_mutect2_bcftools_filter_pass/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz"
        ),
        protected(
            "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz.tbi"
        ),
    threads: 1
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (attempt * 200) + 500,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2i_bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_view/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_bcftools_view",
            default="--with-header --write-index='tbi'",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/view"


"""
Reported on Flamingo on 150 datasets
* time 20s ± 2s
* mem 540mb ± 70mb
"""


rule fair_gatk_mutect2_bcftools_mutect2_stats:
    input:
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz",
        "results/{species}.{build}.{release}.{datatype}/VariantCalling/VCF/{sample}.vcf.gz.tbi",
    output:
        temp(
            "tmp/fair_gatk_mutect2_bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.stats.txt"
        ),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 600 + 100 * attempt,
        runtime=lambda wildcards, attempt: 15 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_mutect2_stats/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        lookup_config(dpath="params/fair_gatk_mutect2_bcftools_stats", default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/stats"

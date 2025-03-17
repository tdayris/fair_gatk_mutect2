"""
## Memory
Requires a job with at most 440.72  Mb,
 on average 328.9 ± 175.33 Mb, 
on Gustave Roussy's HPC Flamingo, on a 91.0  Mb dataset.
## Time
A job took 0:01:46 to proceed,
on average 0:00:39 ± 0:00:38
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
        mem_mb=lambda wildcards, attempt: attempt * 500,
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
        "v5.8.3/bio/bcftools/norm"


"""
## Memory
Requires a job with at most 426.09  Mb,
 on average 316.03 ± 182.31 Mb, 
on Gustave Roussy's HPC Flamingo, on a 35.0  Mb dataset.
## Time
A job took 0:00:36 to proceed,
on average 0:00:16 ± 0:00:13
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
        mem_mb=lambda wildcards, attempt: attempt * 500,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_bcftools_filter_pass/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_bcftools_filter_pass/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra='--include \'FILTER="PASS" || FILTER="."\'',
    wrapper:
        "v5.8.3/bio/bcftools/filter"


"""
## Memory
Requires a job with at most 426.09  Mb,
 on average 316.03 ± 182.31 Mb, 
on Gustave Roussy's HPC Flamingo, on a 35.0  Mb dataset.
## Time
A job took 0:00:36 to proceed,
on average 0:00:16 ± 0:00:13
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
        mem_mb=lambda wildcards, attempt: attempt * 500,
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
        "v5.8.3/bio/bcftools/view"


"""
## Memory
Requires a job with at most 563.22  Mb,
 on average 401.76 ± 235.39 Mb, 
on Gustave Roussy's HPC Flamingo, on a 35.0  Mb dataset.
## Time
A job took 0:00:35 to proceed,
on average 0:00:11 ± 0:00:12
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
        "v5.8.3/bio/bcftools/stats"

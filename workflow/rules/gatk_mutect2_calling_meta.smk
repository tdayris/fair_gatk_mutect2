module gatk_mutect2_calling:
    meta_wrapper:
        f"{snakemake_wrappers_prefix}/meta/bio/gatk_mutect2_calling"
    config:
        config


"""
Reported on Flamingo on 150 datasets
* time 95min ± 49min
* mem  6.5Go ± 3Go
"""


use rule picard_replace_read_groups from gatk_mutect2_calling as fair_gatk_mutect2_picard_reaplace_read_groups with:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        temp(
            "tmp/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4_000 + 6_000,
        runtime=lambda wildcards, attempt: attempt * 45 + 105,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra="--RGPL illumina --RGLB WES_Germline_{species}.{build}.{release}.{datatype} --RGPU {sample} --RGSM {sample}",


"""
Reported on Flamingo on 150 datastets
* time 18min ± 11min (peak at 91min)
* mem 462mb ± 21mb
"""


use rule sambamba_index_picard_bam from gatk_mutect2_calling as fair_gatk_mutect2_sambamba_index_picard_bam with:
    input:
        "tmp/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam.bai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 100 + 500,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_sambamba_index_picard_bam/{species}.{build}.{release}.{datatype}/{sample}_picard_groups.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_sambamba_index_picard_bam/{species}.{build}.{release}.{datatype}/{sample}_picard_groups.tsv"


"""
Reported on Flamingo on 150 datastets
* mem  21Go ± 10Go
* time 5h ± 2h40
"""


use rule mutect2_call from gatk_mutect2_calling as fair_gatk_mutect2_gatk_mutect2_call with:
    input:
        unpack(get_mutect2_call_input),
    output:
        vcf=temp(
            "tmp/fair_gatk_mutect2_gatk_mutect2_call/{species}.{build}.{release}.{datatype}/{sample}.vcf"
        ),
        f1r2=temp(
            "tmp/fair_gatk_mutect2_gatk_mutect2_call/{species}.{build}.{release}.{datatype}/{sample}.f1r2.tar.gz"
        ),
        stats=temp(
            "tmp/fair_gatk_mutect2_gatk_mutect2_call/{species}.{build}.{release}.{datatype}/{sample}.vcf.stats"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 22_000,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_gatk_mutect2_call/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_gatk_mutect2_call/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        use_parallelgc=True,
        use_omp=True,
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_gatk_mutect2_call",
            default=lambda wildcards: (
                ""
                if get_normal_sample(wildcards) is None
                else f" --normal {get_normal_sample(wildcards)}"
            ),
        ),


use rule gatk_get_pileup_summaries from gatk_mutect2_calling as fair_gatk_mutect2_gatk_get_pileup_summaries with:
    input:
        unpack(get_gatk_get_pileup_summaries_input),
    output:
        temp(
            "tmp/fair_gatk_mutect2_gatk_get_pileup_summaries/{species}.{build}.{release}.{datatype}/{sample}.table"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20_000,
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_gatk_get_pileup_summaries{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_gatk_get_pileup_summaries/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_gatk_getpileupsummaries",
            default="",
        ),


use rule gatk_calculate_contamination from gatk_mutect2_calling as fair_gatk_mutect2_gatk_calcultate_contamination with:
    input:
        "tmp/fair_gatk_mutect2_gatk_get_pileup_summaries/{species}.{build}.{release}.{datatype}/{sample}.table",
    output:
        temp(
            "tmp/fair_gatk_mutect2_gatk_calculate_contamination/{species}.{build}.{release}.{datatype}/{sample}.pileups.table"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20_000,
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_gatk_calculate_contamination/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_gatk_calculate_contamination/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_gatk_calculatecontamination",
            default="",
        ),


"""
Reported on Flamingo on 150 datastets
* mem 22Go ± 10Go
* time 1:20 ± 1min
"""


use rule gatk_learn_read_orientation_model from gatk_mutect2_calling as fair_gatk_mutect2_gatk_learn_read_orientation_model with:
    input:
        f1r2="tmp/fair_gatk_mutect2_gatk_mutect2_call/{species}.{build}.{release}.{datatype}/{sample}.f1r2.tar.gz",
    output:
        temp(
            "tmp/fair_gatk_mutect2_gatk_learn_read_orientation_model/{species}.{build}.{release}.{datatype}/{sample}.tar.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8_000 + 22_000,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_gatk_learn_read_orientation_model/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_gatk_learn_read_orientation_model/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_gatk_learnreadorientationmodel",
            default="",
        ),


"""
Reported on Flamingo on ~150 datasets
* time 2h ± 30min
* mem  24Go ± 10Go
"""


use rule filter_mutect_calls from gatk_mutect2_calling as fair_gatk_mutect2_filter_mutect_calls with:
    input:
        unpack(get_filter_mutect_calls_input),
    output:
        vcf=temp(
            "tmp/fair_gatk_mutect2_filter_mutect_calls/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz"
        ),
        vcf_tbi=temp(
            "tmp/fair_gatk_mutect2_filter_mutect_calls/{species}.{build}.{release}.{datatype}/{sample}.vcf.gz.tbi"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10_000 + 24_000,
        runtime=lambda wildcards, attempt: attempt * 60 + 75,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_filtermutectcalls/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_filtermutectcalls/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_gatk_filtermutectcalls",
            default="--create-output-variant-index --min-median-mapping-quality 35 --max-alt-allele-count 3",
        ),

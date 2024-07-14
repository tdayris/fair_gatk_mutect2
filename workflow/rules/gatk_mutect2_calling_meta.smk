module gatk_mutect2_calling:
    meta_wrapper:
        f"{snakemake_wrappers_prefix}/meta/bio/gatk_mutect2_calling"
    config:
        config


"""
## Memory
Requires a job with at most 16386.32  Mb,
 on average 8120.21 ± 5260.7 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:23:24 to proceed,
on average 0:06:13 ± 0:07:33
"""


use rule picard_replace_read_groups from gatk_mutect2_calling as fair_gatk_mutect2_picard_replace_read_groups with:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        temp(
            "tmp/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4_000 + 13_000,
        runtime=lambda wildcards, attempt: attempt * 45 + 85,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra="--RGPL illumina --RGLB WES_Germline_{species}.{build}.{release}.{datatype} --RGPU {sample} --RGSM {sample}",


"""
## Memory
Requires a job with at most 468.36  Mb,
 on average 348.96 ± 183.95 Mb, 
on Gustave Roussy's HPC Flamingo, on a 93.0  Mb dataset.
## Time
A job took 0:05:42 to proceed,
on average 0:01:33 ± 0:01:47
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
## Memory
Requires a job with at most 78070.07  Mb,
 on average 25200.52 ± 23083.55 Mb, 
on Gustave Roussy's HPC Flamingo, on a 92.0  Mb dataset.
## Time
A job took 7:31:03 to proceed,
on average 2:06:11 ± 2:26:14
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
        mem_mb=lambda wildcards, attempt: attempt * 40_000,
        runtime=lambda wildcards, attempt: attempt * 60 * 15,
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
## Memory
Requires a job with at most 69707.14  Mb,
 on average 22480.43 ± 22444.62 Mb, 
on Gustave Roussy's HPC Flamingo, on a 92.0  Mb dataset.

## Time
A job took 0:01:32 to proceed,
on average 0:00:29 ± 0:00:30
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
        mem_mb=lambda wildcards, attempt: attempt * 30_000,
        runtime=lambda wildcards, attempt: attempt * 60,
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
## Memory
Requires a job with at most 56514.84  Mb,
 on average 24296.87 ± 17262.98 Mb, 
on Gustave Roussy's HPC Flamingo, on a 91.0  Mb dataset.
## Time
A job took 0:30:32 to proceed,
on average 0:09:10 ± 0:09:40
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
        mem_mb=lambda wildcards, attempt: attempt * 30_000,
        runtime=lambda wildcards, attempt: attempt * 60
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

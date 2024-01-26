module gatk_mutect2_calling:
    meta_wrapper:
        "v3.3.3/meta/bio/gatk_mutect2_calling"
    config:
        config


use rule picard_replace_read_groups from gatk_mutect2_calling with:
    input:
        "results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
    output:
        temp(
            "tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp",
    log:
        "logs/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra="--RGPL illumina --RGLB WES_Germline_{species}.{build}.{release}.{datatype} --RGPU {sample} --RGSM {sample}",


use rule sambamba_index_picard_bam from gatk_mutect2_calling with:
    input:
        "tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam",
    output:
        temp(
            "tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam.bai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/sambamba/index/{species}.{build}.{release}.{datatype}/{sample}_picard_groups.log",
    benchmark:
        "benchmark/sambamba/index/{species}.{build}.{release}.{datatype}/{sample}_picard_groups.tsv"


use rule mutect2_call from gatk_mutect2_calling with:
    input:
        unpack(get_mutect2_call_input),
    output:
        vcf=temp("tmp/gatk/mutect2/{species}.{build}.{release}.{datatype}/{sample}.vcf"),
        f1r2=temp(
            "tmp/gatk/mutect2/{species}.{build}.{release}.{datatype}/{sample}.f1r2.tar.gz"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/gatk/mutect2/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/gatk/mutect2/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        use_parallelgc=True,
        use_omp=True,
        extra=config.get("params", {}).get("gatk", {}).get("mutect2", ""),


use rule gatk_get_pileup_summaries from gatk_mutect2_calling with:
    input:
        unpack(get_gatk_get_pileup_summaries_input),
    output:
        temp(
            "tmp/gatk/getpileupsummaries/{species}.{build}.{release}.{datatype}/{sample}.table"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/gatk/getpileupsummaries/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/gatk/getpileupsummaries/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {}).get("gatk", {}).get("getpileupsummaries", ""),


use rule gatk_calculate_contamination from gatk_mutect2_calling with:
    input:
        "tmp/gatk/getpileupsummaries/{species}.{build}.{release}.{datatype}/{sample}.table",
    output:
        temp(
            "tmp/gatk/calculatecontamination/{species}.{build}.{release}.{datatype}/{sample}.pileups.table"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/gatk/calculatecontamination/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/gatk/calculatecontamination/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {}).get("gatk", {}).get("calculatecontamination", ""),


use rule gatk_learn_read_orientation_model from gatk_mutect2_calling with:
    input:
        f1r2="tmp/gatk/mutect2/{species}.{build}.{release}.{datatype}/{sample}.f1r2.tar.gz",
    output:
        temp(
            "tmp/gatk/learnreadorientationmodel/{species}.{build}.{release}.{datatype}/{sample}.tar.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/gatk/learnreadorientationmodel/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/gatk/learnreadorientationmodel/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {})
        .get("gatk", {})
        .get("learnreadorientationmodel", ""),


use rule filter_mutect_calls from gatk_mutect2_calling with:
    input:
        unpack(get_filter_mutect_calls_input),
    output:
        vcf="results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz",
        vcf_tbi="results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (20 * 1024),
        runtime=lambda wildcards, attempt: attempt * (60 * 2),
        tmpdir="tmp",
    log:
        "logs/gatk/filtermutectcalls/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/gatk/filtermutectcalls/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {})
        .get("gatk", {})
        .get("filtermutectcalls", "--create-output-variant-index"),

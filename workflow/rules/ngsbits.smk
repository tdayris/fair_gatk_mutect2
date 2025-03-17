rule fair_gatk_mutect2_ngsbits_sampleancestry:
    input:
        lambda wildcards: get_all_vcf_from_genotype(wildcards, samples),
    output:
        report(
            "results/{species}.{build}.{release}.{datatype}/VariantCalling/QC/Ancestry.tsv",
            caption="../report/ngsbits_ancestry.rst",
            category="Quality Control",
            subcategory="General",
            labels={
                "report": "tsv",
                "step": "Mapping",
                "organism": "{species}.{build}.{release}.{datatype}",
            },
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_ngsbits_sampleancestry/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_ngsbits_sampleancestry/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lambda wildcards: (
            " -build hg38 "
            if str(wildcards.build).lower() == "grch38"
            else " -build hg19 "
        ),
    wrapper:
        "v5.8.3/bio/ngsbits/sampleancestry"


rule fair_gatk_mutect2_ngsbits_samplesimilarity:
    input:
        ref=lambda wildcards: select_fai(wildcards),
        regions=lambda wildcards: get_intervals(wildcards),
        samples=get_all_vcf_from_genotype(wildcards, samples),
    output:
        report(
            "results/{species}.{build}.{release}.{datatype}/VariantCalling/QC/Similarity.tsv",
            caption="../report/ngsbits_similarity.rst",
            category="Quality Control",
            subcategory="General",
            labels={
                "report": "tsv",
                "step": "Variant Calling",
                "organism": "{species}.{build}.{release}.{datatype}",
            },
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_gatk_mutect2_ngsbits_samplesimilarity/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_ngsbits_samplesimilarity/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lambda wildcards: (
            " -build hg38 "
            if str(wildcards.build).lower() == "grch38"
            else " -build hg19 "
        ),
    wrapper:
        "v5.8.3/bio/ngsbits/samplesimilarity"

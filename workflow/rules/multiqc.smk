rule fair_gatk_mutect_germline_multiqc_report:
    input:
        picard_qc=collect(
            "tmp/fair_bowtie2_mapping/picard_create_multiple_metrics/{sample.species}.{sample.build}.{sample.release}.dna/stats/{sample.sample_id}{ext}",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            ext=[
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".base_distribution_by_cycle_metrics",
                ".base_distribution_by_cycle.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
            ],
        ),
        fastp_pair_ended=collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastp_single_ended=collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastqc_pair_ended=collect(
            "results/QC/report_pe/{sample.sample_id}.{stream}_fastqc.zip",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            stream=stream_list,
        ),
        fastqc_single_ended=collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bowtie2=collect(
            "logs/fair_bowtie2_mapping/bowtie2_alignment/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.log",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        samtools=collect(
            "tmp/fair_bowtie2_mapping/samtools_stats/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bcftools_log=collect(
            "logs/fair_gatk_mutect_germline/bcftools_mutect2_stats/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.log",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bcftools_stats=collect(
            "tmp/fair_gatk_mutect_germline/bcftools_mutect2_stats/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.stats.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
    output:
        report(
            "results/{species}.{build}.{release}.dna/QC/MultiQC_GatkGermlineCalling.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "GatkGermlineCalling",
                "organirm": "{species}.{build}.{release}",
            },
        ),
        "results/{species}.{build}.{release}.dna/QC/MultiQC_GatkGermlineCalling_data.zip",
    params:
        extra=lookup(dpath="params/multiqc", within=config),
        use_input_files_only=True,
    log:
        "logs/fair_gatk_mutect_germline/multiqc_report/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/multiqc_report/{species}.{build}.{release}.tsv"
    wrapper:
        "v3.3.6/bio/multiqc"

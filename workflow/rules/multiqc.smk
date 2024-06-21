"""
Reported on Flamingo
* time 10s ± 2s
* mem  356mb ± 1mb
"""


rule fair_gatk_mutect2_multiqc_config:
    input:
        "tmp/fair_fastqc_multiqc_bigr_logo.png",
    output:
        temp("tmp/fair_gatk_mutect2_multiqc_config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 100 + 370,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    localrule: True
    log:
        "logs/fair_gatk_mutect2_multiqc_config.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_multiqc_config.tsv"
    params:
        extra=lambda wildcards, input: {
            "title": "Variant Calling quality control report",
            "subtitle": "Produced on raw fastq recieved from sequencer",
            "intro_text": (
                "This pipeline building this report has "
                "no information about sequencing protocol, "
                "nor wet-lab experimental design."
            ),
            "report_comment": (
                "This report was generated using: "
                "https://github.com/tdayris/fair_gatk_mutect2"
            ),
            "show_analysis_paths": False,
            "show_analysis_time": False,
            "custom_logo": input[0],
            "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
            "custom_logo_title": "Bioinformatics Platform @ Gustave Roussy",
            "report_header_info": [
                {"Contact E-mail": "bigr@gustaveroussy.fr"},
                {"Application type": "Short-gapped reads"},
                {"Project Type": "Whole Exome calling"},
            ],
            "software_versions": {
                "Quality controls": {
                    "fastqc": "1.12.1",
                    "fastq_screen": "0.15.3",
                    "bowtie2": "1.3.1",
                    "multiqc": "1.20.0",
                },
                "Mapping": {
                    "bowtie2": "2.5.3",
                    "sambamba": "1.0",
                    "samtools": "1.19.2",
                    "picard": "3.1.1",
                    "rseqc": "5.0.3",
                    "fastp": "0.23.4",
                    "ngsderive": "3.3.2",
                    "goleft": "0.2.4",
                },
                "Variant Calling": {
                    "picard": "3.1.1",
                    "gatk": "4.5.0.0",
                    "bcftools": "1.19",
                    "snpeff": "5.2",
                    "snpsift": "5.2",
                },
                "Genome Indexing": {
                    "agat": "1.3.1",
                    "picard": "3.1.1",
                    "gffread": "0.12.7",
                    "samtools": "1.0",
                },
                "Pipeline": {
                    "snakemake": "8.5.3",
                    "fair_fastqc_multiqc": "2.1.2",
                    "fair_genome_indexer": "3.2.2",
                    "fair_bowtie2_mapping": "3.2.0",
                    "fair_gatk_mutect2": "1.2.0",
                },
            },
            "disable_version_detection": True,
            "run_modules": [
                "fastqc",
                "fastq_screen",
                "fastp",
                "bowtie2",
                "samtools",
                "picard",
                "rseqc",
                "ngsderive",
                "goleft_indexcov",
                "gatk",
                "bcftools",
                "snpeff",
            ],
            "report_section_order": {
                "fastq_screen": {"order": 1000},
                "ngsderive": {"order": 950},
                "fastqc": {"order": 900},
                "fastp": {"order": 890},
                "bowtie2": {"order": 880},
                "picard": {"order": 870},
                "samtools": {"order": 860},
                "rseqc": {"order": 850},
                "goleft_indexcov": {"order": 840},
                "gatk": {"order": 835},
                "bcftools": {"order": 830},
                "snpeff": {"order": 820},
                "software_versions": {"order": -1000},
            },
        },
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_gatk_mutect2_multiqc_config.py"


rule fair_gatk_mutect2_multiqc_report:
    input:
        config="tmp/fair_gatk_mutect2_multiqc_config.yaml",
        picard_qc=collect(
            "tmp/fair_bowtie2_mapping_picard_create_multiple_metrics/{sample.species}.{sample.build}.{sample.release}.dna/stats/{sample.sample_id}{ext}",
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
            "tmp/fair_bowtie2_mapping_fastp_trimming_pair_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastp_single_ended=collect(
            "tmp/fair_bowtie2_mapping_fastp_trimming_single_ended/{sample.sample_id}.fastp.json",
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
            stream=stream_tuple,
        ),
        fastqc_single_ended=collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bowtie2=collect(
            "logs/fair_bowtie2_mapping_bowtie2_alignment/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.log",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        samtools=collect(
            "tmp/fair_bowtie2_mapping_samtools_stats/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_infer_experiment=collect(
            "tmp/fair_bowtie2_mapping_rseqc_infer_experiment/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.infer_experiment.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_bamstat=collect(
            "tmp/fair_bowtie2_mapping_rseqc_bamstat/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.bamstat.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_read_gc=collect(
            "tmp/fair_bowtie2_mapping_rseqc_read_gc/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.GC.xls",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_read_distribution=collect(
            "tmp/fair_bowtie2_mapping_rseqc_read_distribution/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_inner_distance=collect(
            "tmp/fair_bowtie2_mapping_rseqc_inner_distance/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.inner_distance_freq.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        goleft_indexcov_ped=collect(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{sample.species}.{sample.release}.{sample.build}.dna/{sample.sample_id}-indexcov.ped",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        goleft_indexcov_roc=collect(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{sample.species}.{sample.release}.{sample.build}.dna/{sample.sample_id}-indexcov.roc",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        ngsderive_readlen=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_readlen/{sample.species}.{sample.build}/{sample.release}.dna/{sample.sample_id}.readlen.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        ngsderive_instrument=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_instrument/{sample.species}.{sample.build}/{sample.release}.dna/{sample.sample_id}.instrument.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        ngsderive_encoding=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_encoding/{sample.species}.{sample.build}/{sample.release}.dna/{sample.sample_id}.encoding.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bcftools_log=collect(
            "logs/fair_gatk_mutect2_bcftools_mutect2_stats/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.log",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bcftools_stats=collect(
            "tmp/fair_gatk_mutect2_bcftools_mutect2_stats/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.stats.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        snpeff_stats=collect(
            "tmp/fair_gatk_mutect2_snpeff_annotate/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.csv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        variant_eval=collect(
            "tmp/fair_gatk_mutect2_gatk_germline_varianteval/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.grp",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
    output:
        report(
            "results/{species}.{build}.{release}.dna/QC/MultiQC_GatkCalling.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "GatkCalling",
                "organirm": "{species}.{build}.{release}",
            },
        ),
        "results/{species}.{build}.{release}.dna/QC/MultiQC_GatkCalling_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    params:
        extra=lookup_config(
            dpath="params/fair_gatk_mutect2_multiqc",
            default="--verbose --no-megaqc-upload --no-ansi --force",
        ),
        use_input_files_only=True,
    log:
        "logs/fair_gatk_mutect2_multiqc_report/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_gatk_mutect2_multiqc_report/{species}.{build}.{release}.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/multiqc"

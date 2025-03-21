module fair_bowtie2_mapping:
    snakefile:
        config.get(
            "fair_bowtie2_mapping",
            github(
                "tdayris/fair_bowtie2_mapping",
                path="workflow/Snakefile",
                tag="4.1.3",
            ),
        )
    config:
        {
            **config,
            "load_fair_genome_indexer": False,
            "load_fair_fastqc_multiqc": False,
        }


use rule * from fair_bowtie2_mapping

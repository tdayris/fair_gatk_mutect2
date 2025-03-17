module fair_star_mapping:
    snakefile:
        config.get(
            "fair_star_mapping",
            github(
                "tdayris/fair_star_mapping",
                path="workflow/Snakefile",
                tag="1.1.0",
            ),
        )
    config:
        {
            **config,
            "load_fair_genome_indexer": False,
            "load_fair_fastqc_multiqc": False,
            "load_bowtie2_mapping": True,
        }


use rule * from fair_star_mapping

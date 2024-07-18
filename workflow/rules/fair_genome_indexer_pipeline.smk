module fair_genome_indexer:
    snakefile:
        config.get(
            "fair_genome_indexer_snakefile",
            github(
                "tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="3.8.1"
            ),
        )
    config:
        config


use rule * from fair_genome_indexer

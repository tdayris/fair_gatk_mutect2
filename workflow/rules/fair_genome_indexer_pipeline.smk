module fair_genome_indexer:
    snakefile:
        github("tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="3.1.4")
    config:
        config


use rule * from fair_genome_indexer

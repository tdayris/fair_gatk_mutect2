module fair_fastqc_multiqc:
    snakefile:
        github("tdayris/fair_fastqc_multiqc", path="workflow/Snakefile", tag="1.0.4")
    config:
        {
            "samples": config.get("samples", "config/samples.csv"),
            "genomes": config.get("genomes", "config/genomes.csv"),
            "params": {
                "fastqc": config.get("params", {}).get("fastqc", ""),
                "multiqc": config.get("params", {}).get("multiqc", "--zip-data-dir"),
            },
        }


use rule * from fair_fastqc_multiqc

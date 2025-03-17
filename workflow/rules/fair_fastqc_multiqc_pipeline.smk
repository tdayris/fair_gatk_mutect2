module fair_fastqc_multiqc:
    snakefile:
        config.get(
            "fair_fastqc_multiqc",
            github(
                "tdayris/fair_fastqc_multiqc",
                path="workflow/Snakefile",
                tag="2.3.5",
            ),
        )
    config:
        config


use rule * from fair_fastqc_multiqc

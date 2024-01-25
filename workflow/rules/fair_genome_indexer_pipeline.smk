module fair_genome_indexer:
    snakefile:
        github("tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="3.0.0")
    config:
        {
            "genomes": config.get("genomes", "config/genomes.csv"),
            "params": {
                "pyfaidx": config.get("params", {}).get(
                    "pyfaidx", {"dna": '--regex "^[[0-9]+|X|Y|MT]"', "cdna": ""}
                ),
                "samtools": config.get("params", {}).get("samtools", {"faidx": ""}),
                "picard": config.get("params", {}).get(
                    "picard", {"createsequencedictionary": ""}
                ),
                "tabix": config.get("params", {}).get("tabix", "-p vcf"),
                "bedtools": config.get("params", {}).get(
                    "bedtools", {"filter_non_canonical_chrom": ""}
                ),
                "gffread": config.get("params", {}).get("gffread", ""),
                "agat": config.get("params", {}).get(
                    "agat",
                    {
                        "gff2gtf": "",
                        "filter_features": "",
                        "agat_convert_sp_gff2tsv": "",
                        "select_feature_by_attribute_value": "--attribute 'transcript_support_level' --value '\"NA\"' --test '='",
                    },
                ),
            },
        }


use rule * from fair_genome_indexer

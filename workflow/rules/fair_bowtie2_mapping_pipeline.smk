module fair_bowtie2_mapping:
    snakefile: github("tdayris/fair_bowtie2_mapping", path="workflow/Snakefile", tag="3.0.0")
    config: 
        {
            "samples": config.get("samples", "config/samples.csv"),
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
                "fastp": config.get("params", {}).get("fastp", {"adapters": "", "extra": ""}),
                "fastqc": config.get("params", {}).get("fastqc", ""),
                "bowtie2": config.get("params", {}).get("bowtie2", {"build": "", "align": ""}),
                "sambamba": {
                    "view": config.get("params", {})
                    .get("sambamba", {})
                    .get("view", "--format 'bam'"),
                    "markdup": config.get("params", {})
                    .get("sambamba", {})
                    .get("markdup", "--overflow-list-size=500000"),
                },
                "picard": config.get("params", {}).get("picard", {"metrics": "", "createsequencedictionary": ""}),
                "samtools": config.get("params", {}).get("samtools", {"faidx": "", "stats": ""}),
                "multiqc": config.get("params", {}).get(
                    "multiqc", 
                    "--module deeptools --module macs2 --module picard --module fastqc "
                    "--module fastp --module samtools --module bowtie2 --module sambamba "
                    "--zip-data-dir --verbose --no-megaqc-upload --no-ansi --force"
                ),
            },
            "genomes": config.get("genomes", "config/genomes.csv"),
            "load_fair_genome_indexer": False,
            "load_fair_fastqc_multiqc": False,
        }

use rule * from fair_bowtie2_mapping
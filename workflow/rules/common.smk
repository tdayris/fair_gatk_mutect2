import csv
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from typing import Any

snakemake.utils.min_version("7.29.0")

# containerized: "docker://snakemake/snakemake:v7.32.4"
# containerized: "docker://mambaorg/micromamba:git-8440cec-jammy-cuda-12.2.0"
# containerized: "docker://condaforge/mambaforge:23.3.1-1"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check samples properties table
sample_table_path: str = config.get("samples", "config/samples.csv")
with open(sample_table_path, "r") as sample_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(sample_table_stream.readline())
    sample_table_stream.seek(0)

samples: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=sample_table_path,
    sep=dialect.delimiter,
    header=0,
    index_col=None,
    comment="#",
    dtype=str,
)
samples = samples.where(samples.notnull(), None)
snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")

# This is here for compatibility with
genome_table_path: str = config.get("genomes")
if genome_table_path:
    with open(genome_table_path, "r") as genome_table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(genome_table_stream.readline())
        genome_table_stream.seek(0)

    genomes: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=genome_table_path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )
    genomes = genomes.where(genomes.notnull(), None)
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflows.rst"


release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
datatype_list: list[str] = ["dna", "cdna", "transcripts"]
stream_list: list[str] = ["1", "2"]


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),
    stream=r"|".join(stream_list),


def get_reference_genome_data(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return genome information for a given set of {species, build, release} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str | None]):
    Genome information
    """
    result: str | None = genomes.loc[
        (genomes["species"] == str(wildcards.species))
        & (genomes["build"] == str(wildcards.build))
        & (genomes["release"] == str(wildcards.release))
    ]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_sample_information(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return sample information for a given {sample} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their input files

    Return (dict[str, str | None]):
    Sample information
    """
    result: str | None = samples.loc[(samples["sample_id"] == str(wildcards.sample))]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_mutect2_call_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
):
    """
    Return best input files list, according to
    Mutect2-snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionary of input files
    """
    genome_data: dict[str, str | None] = get_reference_genome_data(wildcards, genomes)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    mutect2_call_input: dict[str, str] = {
        "fasta": genome_data.get(
            "dna_fasta",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
        ),
        "fasta_fai": genome_data.get(
            "dna_fai",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
        ),
        "fasta_dict": genome_data.get(
            "dna_dict",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.dict",
        ),
        "map": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam",
        "map_bai": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
    }

    capture_kit: str | None = genome_data.get("capture_kit")
    if capture_kit:
        mutect2_call_input["intervals"] = capture_kit

    pon: str | None = genome_data.get("PoN")
    if pon:
        mutect2_call_input["pon"] = pon

    af_only: str | None = genome_data.get("af_only_vcf")
    af_only_tbi: str | None = genome_data.get("af_only_tbi")
    if af_only and af_only_tbi:
        mutect2_call_input["germline"] = af_only
        mutect2_call_input["germline_tbi"] = af_only_tbi

    return mutect2_call_input


def get_gatk_get_pileup_summaries_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
):
    """
    Return best input files list, according to
    GATK-GetPileupSummaries snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionary of input files
    """
    genome_data: dict[str, str | None] = get_reference_genome_data(wildcards, genomes)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    gatk_get_pileup_summaries_input: dict[str, str] = {
        "bam": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam",
        "bam_bai": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
    }

    capture_kit: str | None = genome_data.get("capture_kit")
    if capture_kit:
        gatk_get_pileup_summaries_input["intervals"] = capture_kit

    af_only: str | None = genome_data.get("af_only_vcf")
    af_only_tbi: str | None = genome_data.get("af_only_tbi")
    if af_only and af_only_tbi:
        gatk_get_pileup_summaries_input["variants"] = af_only
        gatk_get_pileup_summaries_input["variants_tbi"] = af_only_tbi

    return gatk_get_pileup_summaries_input


def get_filter_mutect_calls_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
):
    """
    Return best input files list, according to
    GATK-FilterMutectCall snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionary of input files
    """
    genome_data: dict[str, str | None] = get_reference_genome_data(wildcards, genomes)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    filter_mutect_calls_input: dict[str, str] = {
        "ref": genome_data.get(
            "dna_fasta",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
        ),
        "ref_fai": genome_data.get(
            "dna_fai",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
        ),
        "ref_dict": genome_data.get(
            "dna_dict",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.dict",
        ),
        "aln": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam",
        "aln_idx": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
        "vcf": f"tmp/gatk/mutect2/{species}.{build}.{release}.{datatype}/{sample}.vcf",
        "contamination": f"tmp/gatk/calculatecontamination/{species}.{build}.{release}.{datatype}/{sample}.pileups.table",
        "f1r2": f"tmp/gatk/learnreadorientationmodel/{species}.{build}.{release}.{datatype}/{sample}.tar.gz",
    }

    capture_kit: str | None = genome_data.get("capture_kit")
    if capture_kit:
        filter_mutect_calls_input["intervals"] = capture_kit

    return filter_mutect_calls_input


def get_gatk_germline_varianteval_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
):
    """
    Return best input files list, according to
    GATK-FilterMutectCall snakemake-wrappers' requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str]):
    Dictionary of input files
    """
    genome_data: dict[str, str | None] = get_reference_genome_data(wildcards, genomes)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    return {
        "vcf": f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz",
        "vcf_tbi": f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz.tbi",
        "bam": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam",
        "bai": f"tmp/picard/add_or_replace_groupe/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
        "ref": genome_data.get(
            "dna_fasta",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
        ),
        "dict": genome_data.get(
            "dna_dict",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.dict",
        ),
        "fai": genome_data.get(
            "dna_fai",
            f"reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
        ),
        "known": genome_data.get(
            "dbsnp", f"reference/variants/{species}.{build}.{release}.all.vcf.gz",
        ),
        "known_tbi": genome_data.get(
            "dbsnp_tbi",
            f"reference/variants/{species}.{build}.{release}.all.vcf.gz.tbi",
        ),
    }


def get_fair_gatk_mutect_germline_multiqc_report_input(
    wildcards: snakemake.io.Wildcards,
    sample: pandas.DataFrame = samples,
) -> dict[str, str | list[str]]:
    """
    Return expected input files for Multiqc peak calling report,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, Union[str, list[str]]]):
    Input files dict, as required by MultiQC's snakemake-wrapper
    """
    results: dict[str, list[str]] = {
        "picard_qc": [],
        "fastp": [],
        "fastqc": [],
        "bowtie2": [],
        "samtools": [],
        "bcftools": [],
        "varianteval": [],
    }
    datatype: str = "dna"
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        results["picard_qc"] += multiext(
            f"tmp/picard/{species}.{build}.{release}.{datatype}/stats/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        )
        sample_data: dict[str, str | None] = get_sample_information(
            snakemake.io.Wildcards(fromdict={"sample": sample}), samples
        )
        if sample_data.get("downstream_file"):
            results["fastp"].append(f"tmp/fastp/report_pe/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.1_fastqc.zip")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.2_fastqc.zip")
        else:
            results["fastp"].append(f"tmp/fastp/report_se/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}_fastqc.zip")

        results["bowtie2"].append(
            f"logs/bowtie2/align/{species}.{build}.{release}.{datatype}/{sample}.log"
        )

        results["samtools"].append(
            f"tmp/samtools/{species}.{build}.{release}.{datatype}/{sample}.txt"
        )

        results["bcftools"].append(
            f"tmp/bcftools/stats/mutect2/germline/{species}.{build}.{release}.{datatype}/{sample}.stats.txt"
        )

        results["varianteval"].append(
            f"tmp/gatk/varianteval/{species}.{build}.{release}.{datatype}/{sample}.grp"
        )

    return results


def get_gatk_mutect_germline_targets(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
) -> dict[str, list[str]]:
    """
    Return expected output files for this pipeline,
    according to user-input, and snakemake requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, list[str]]):
    Output files dict
    """
    results: dict[str, list[str]] = {
        "multiqc": [
            "results/QC/MultiQC_FastQC.html",
            "results/QC/MultiQC_Mapping.html",
            "results/QC/MultiQC_GatkGermlineCalling.html",
        ],
        "vcf": [],
        "vcf_tbi": [],
    }
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    datatype: str = "dna"

    for sample, species, build, release in sample_iterator:
        results[
            "vcf"
        ] = f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz"
        results[
            "vc_tbi"
        ] = f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz.tbi"

    return results

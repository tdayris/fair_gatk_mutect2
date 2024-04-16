import csv
import os
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from typing import Any, NamedTuple

snakemake.utils.min_version("8.4.0")


container: "docker://snakemake/snakemake:v8.5.3"


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
tmp: str = f"{os.getcwd()}/tmp"
snakemake_wrappers_prefix: str = "v3.7.0"


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),
    stream=r"|".join(stream_list),


def lookup_config(
    dpath: str, default: str | None = None, config: dict[str, Any] = config
) -> str:
    """
    Run lookup function with default parameters in order to search a key in configuration and return a default value
    """
    value: str | None = default

    try:
        value = lookup(dpath=dpath, within=config)
    except LookupError:
        value = default
    except WorkflowError:
        value = default

    return value


def lookup_genomes(
    wildcards: snakemake.io.Wildcards,
    key: str,
    default: str | list[str] | None = None,
    genomes: pandas.DataFrame = genomes,
) -> str:
    """
    Run lookup function with default parameters in order to search user-provided sequence/annotation files
    """
    query: str = (
        "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}'".format(
            wildcards=wildcards
        )
    )
    return getattr(lookup(query=query, within=genomes), key, default)


def get_dna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fasta", default=default, genomes=genomes)


def get_cdna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fasta", default=default, genomes=genomes)


def get_transcripts_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fasta", default=default, genomes=genomes
    )


def select_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fasta(wildcards),
            "cdna": get_cdna_fasta(wildcards),
            "transcripts": get_transcripts_fasta(wildcards),
        },
    )


def get_dna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences index
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fai", default=default, genomes=genomes)


def get_cdna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences index
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fai", default=default, genomes=genomes)


def get_transcripts_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences index
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fai", default=default, genomes=genomes
    )


def get_dna_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.dict".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_dict", default=default, genomes=genomes)


def get_cdna_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.dict".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_dict", default=default, genomes=genomes)


def get_transcripts_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.dict".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_dict", default=default, genomes=genomes
    )


def select_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcards and return the right fasta dictionary
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_dict(wildcards),
            "cdna": get_cdna_dict(wildcards),
            "transcripts": get_transcripts_dict(wildcards),
        },
    )


def select_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta index file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fai(wildcards),
            "cdna": get_cdna_fai(wildcards),
            "transcripts": get_transcripts_fai(wildcards),
        },
    )


def get_gtf(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation
    """
    default: str = (
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}.gtf".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="gtf", default=default, genomes=genomes)


def get_intervals(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | None:
    """
    Return path to capturekit file
    """
    return lookup_genomes(wildcards, key="capture_kit", default=None, genomes=genomes)


def get_known_variants(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | None:
    """
    Return path to known variants (AF only VCF)
    """
    return lookup_genomes(wildcards, key="af_only_vcf", default=None, genomes=genomes)


def get_known_variants_tbi(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | None:
    """
    Return path to known variants index file (AF only VCF)
    """
    return lookup_genomes(wildcards, key="af_only_tbi", default=None, genomes=genomes)


def get_pon(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | None:
    """
    Return panel of normals
    """
    return lookup_genomes(wildcards, key="pon", default=None, genomes=genomes)


def get_normal_sample(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> str | None:
    """
    Return corresponding Normal sample (if any)
    """
    query: str = (
        "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}' & sample_id == '{wildcards.sample}'".format(
            wildcards=wildcards
        )
    )
    sample_query: NamedTuple = lookup(query=query, within=samples)
    return getattr(sample_query, "normal_sample_id", None)


def get_normal_bam(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> str | None:
    """
    Return corresponding Normal bam file (if any)
    """
    normal_id: str | None = get_normal_sample(wildcards, samples)
    if normal_id:
        return (
            "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{normal_id.sample_id}.bam".format(
                wildcards=wildcards, normal_id=normal_id
            ),
        )


def get_normal_bai(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> str | None:
    """
    Return corresponding Normal bam index file (if any)
    """
    normal_id: str | None = get_normal_sample(wildcards, sample)
    if normal_id:
        return (
            "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{normal_id.sample_id}.bam.bai".format(
                wildcards=wildcards, normal_id=normal_id
            ),
        )


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
    mutect2_call_input: dict[str, str] = {
        "fasta": get_dna_fasta(wildcards),
        "fasta_fai": get_dna_fai(wildcards),
        "fasta_dict": get_dna_dict(wildcards),
        "map": "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.bam".format(
            wildcards=wildcards
        ),
        "map_bai": "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.bam.bai".format(
            wildcards=wildcards
        ),
    }

    intervals: str | None = get_intervals(wildcards)
    if intervals:
        mutect2_call_input["intervals"] = intervals

    af_only_vcf: str | None = get_known_variants(wildcards)
    af_only_tbi: str | None = get_known_variants_tbi(wildcards)
    if af_only_vcf and af_only_tbi:
        mutect2_call_input["germline"] = af_only_vcf
        mutect2_call_input["germline_tbi"] = af_only_tbi

    pon: str | None = get_pon(wildcards)
    if pon:
        mutect2_call_input["pon"] = pon

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
    gatk_get_pileup_summaries_input: dict[str, str] = {
        "bam": "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.bam".format(
            wildcards=wildcards
        ),
        "bam_bai": "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.bam.bai".format(
            wildcards=wildcards
        ),
    }

    intervals: str | None = get_intervals(wildcards)
    if intervals:
        gatk_get_pileup_summaries_input["intervals"] = intervals

    af_only: str | None = get_known_variants(wildcards)
    af_only_tbi: str | None = get_known_variants_tbi(wildcards)
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
    filter_mutect_calls_input: dict[str, str] = {
        "ref": get_dna_fasta(wildcards),
        "ref_fai": get_dna_fai(wildcards),
        "ref_dict": get_dna_dict(wildcards),
        "aln": "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.bam".format(
            wildcards=wildcards
        ),
        "aln_idx": "tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.bam.bai".format(
            wildcards=wildcards
        ),
        "vcf": "tmp/fair_gatk_mutect2/gatk_mutect2_call/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.vcf".format(
            wildcards=wildcards
        ),
        "f1r2": "tmp/fair_gatk_mutect2/gatk_learn_read_orientation_model/{wildcards.species}.{wildcards.build}.{wildcards.release}.{wildcards.datatype}/{wildcards.sample}.tar.gz".format(
            wildcards=wildcards
        ),
    }

    intervals: str | None = get_intervals(wildcards)
    if intervals:
        gatk_get_pileup_summaries_input["intervals"] = intervals

    af_only: str | None = get_known_variants(wildcards)
    af_only_tbi: str | None = get_known_variants_tbi(wildcards)
    if af_only and af_only_tbi:
        filter_mutect_calls_input["contamination"] = (
            f"tmp/fair_gatk_mutect2/gatk_calcultate_contamination/{species}.{build}.{release}.{datatype}/{sample}.pileups.table"
        )

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
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"
    sample: str = str(wildcards.sample)

    gatk_germline_varianteval_input: dict[str, str] = {
        "vcf": f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz",
        "vcf_tbi": f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Raw/{sample}.vcf.gz.tbi",
        "bam": f"tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam",
        "bai": f"tmp/fair_gatk_mutect2/picard_reaplace_read_groups/{species}.{build}.{release}.{datatype}/{sample}.bam.bai",
        "ref": get_dna_fasta(wildcards),
        "dict": get_dna_dict(wildcards),
        "fai": get_dna_fai(wildcards),
    }

    af_only: str | None = get_known_variants(wildcards)
    af_only_tbi: str | None = get_known_variants_tbi(wildcards)
    if af_only and af_only_tbi:
        gatk_germline_varianteval_input["known"] = af_only
        gatk_germline_varianteval_input["known_tbi"] = af_only_tbi

    return gatk_germline_varianteval_input


def get_gatk_mutect2_targets(
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
        results["multiqc"].append(
            f"results/{species}.{build}.{release}.dna/QC/MultiQC_Mapping.html"
        )
        results["multiqc"].append(
            f"results/{species}.{build}.{release}.dna/QC/MultiQC_GatkGermlineCalling.html"
        )
        results["vcf"].append(
            f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz"
        )
        results["vcf_tbi"].append(
            f"results/{species}.{build}.{release}.{datatype}/VariantCalling/Germline/VCF/{sample}.vcf.gz.tbi"
        )

    results["multiqc"] = list(set(results["multiqc"]))

    return results

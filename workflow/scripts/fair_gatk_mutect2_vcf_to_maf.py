# coding: utf-8

import bz2
import pandas

from collections import defaultdict

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

InfoType = str | int | bool | None | list[str]
AnnType = dict[str, str | None | list[str]]

headers_fmt: tuple[str] = (
    "Genotypes",
    "Allele_Depth",
    "Allele_Fraction",
    "Approximate_Read_Depth",
    "F1R2",
    "F2R1",
    "Fragment_Supporting_Each_Alleles",
    "SB",
)

headers_NMD: tuple[str] = (
    "Gene_Name",
    "Gene_ID",
    "Number_of_transcripts_in_gene",
    "Percent_of_transcripts_affected",
)

headers_ANN: tuple[str] = (
    "Allele",
    "Annotation",
    "Annotation_Impact",
    "Gene_Name",
    "Gene_ID",
    "Feature_Type",
    "Feature_ID",
    "Transcript_BioType",
    "Rank",
    "HGVS.c",
    "HGVS.p",
    "cDNA.pos",
    "CDS.pos",
    "AA.pos",
    "AA.length",
    "Distance",
    "ERRORS_WARNINGS_INFO",
)


def parse_ann(annotation: str, headers: tuple[str] = headers_ANN) -> AnnType:
    """
    Parse ANN info field in VCF, from SnpEff annotation tool
    """
    result = defaultdict(list)
    ann_transcripts: list[str | None] = annotation.split(",")
    for transcript in ann_transcripts:
        values: list[str | None] = transcript.split("|")
        for key, value in zip(headers, values):
            result[key].append(value)

    for key, val in result.items():
        if isinstance(val, list):
            result[key] = ",".join(list(set(val)))

    return result


def info_to_dict(
    info_fields: str, headers_ANN: tuple[str] = headers_ANN
) -> dict[str, InfoType]:
    info_dict: dict[str, InfoType] = {}
    info_fields_list: list[str] = info_fields.split(";")
    for info_field in info_fields_list:
        if "=" in info_field:
            info_key_val: list[str] = info_field.split("=")
            key: str = info_key_val[0]
            val: str = "=".join(info_key_val[1:])

            if key == "ANN":
                info_dict.update(annotation=parse_ann(val), headers=headers_ANN)

            else:
                info_dict[key] = value
        else:
            info_dict[key] = True

    return info_dict


def format_to_dict(
    format_fields: str,
    ncbi_build: str,
    tumor_id: str,
    normal_id: str | None = None,
    headers: tuple[str] = headers_fmt,
) -> dict[str, str | float]:
    """
    Parse format from default GATK keys
    """
    format_list: list[str] = format_fields.split(":")
    return dict(zip(headers, format_list))


def line_to_maf(
    line: str,
    tumod_id: str,
    headers_ANN: tuple[str] = headers_ANN,
    headers_fmt: tuple[str] = headers_fmt,
    normal_id: str | None = None,
) -> dict[str, InfoType]:
    line_fields: list[str] = line.split("\t")
    info_dict: dict[str, InfoType] = info_to_dict(
        info_fields=line_fields[7], headers=headers_ANN
    )
    fmt_dict: dict[str, str | float] = format_to_dict(
        format_fields=line_fields[-1],
        tumod_id=tumod_id,
        normal_id=normal_id,
        headers=headers_fmt,
    )

    result: dict[str, InfoType] = {
        "Hugo_Symbol": info_dict["Gene_Name"],
        "Entrez_Gene_Id": info_dict["Gene_ID"],
        "Center": "Gustave_Roussy",
        "NCBI_Build": ncbi_build,
        "Chromosome": line_fields[0],
        "Start_Position": int(line_fields[1]),
        "End_Position": int(line_fields[1]),
        "Variant_Classification": info_dict["Annotation"],
        "Variant_Type": info_dict["VARTYPE"],
        "Reference_Allele": line_fields[3],
        "Tumor_Seq_Allele1": (
            line_fields[3] if "0" in fmt_dict[tumod_id]["Genotypes"] else line_fields[4]
        ),
        "Tumor_Seq_Allele2": line_fields[4],
        "Tumor_Sample_Barcode": tumor_id,
        "HGVSc": info_dict["HGVS.c"],
        "HGVSp": info_dict["HGVS.p"],
        "Transcript_ID": info_dict["Feature_ID"],
        "t_depth": fmt_dict[tumor_id]["Approximate_Read_Depth"],
        "t_ref_count": fmt_dict[tumor_id]["Fragment_Supporting_Each_Alleles"].split(
            ","
        )[0],
        "t_alt_count": fmt_dict[tumor_id]["Fragment_Supporting_Each_Alleles"].split(
            ","
        )[-1],
        "Gene": info_dict["Gene_Name"],
        "Feature": info_dict["Feature_ID"],
        "Feature_type": info_dict["Feature_Type"],
        "cDNA_position": info_dict["cDNA.pos"],
        "CDS_position": info_dict["CDS.pos"],
        "Protein_position": info_dict["AA.pos"],
        "DISTANCE": info_dict["Distance"],
        "SYMBOL": info_dict["Gene_Name"],
        "SYMBOL_SOURCE": "Ensembl",
        "HGNC_ID": info_dict["Gene_Name"],
        "BIOTYPE": info_dict["Transcript_BioType"],
        "FILTER": line_fields[6],
        "vcf_info": line_fields[7],
        "vcf_format": " ".join(line_fields[8:]),
        "vcf_tumor_gt": fmt_dict[tumor_id]["Genotypes"],
        "IMPACT": info_dict["Annotation_Impact"],
    }

    if normal_id is not None:
        result["Matched_Norm_Sample_Barcode"] = normal_id
        result["Match_Norm_Seq_Allele1"] = line_fields[3]
        result["Match_Norm_Seq_Allele2"] = (
            line_fields[4]
            if "1" in fmt_dict[normal_id]["Genotypes"]
            else line_fields[3]
        )
        result["n_depth"] = fmt_dict[normal_id]["Approximate_Read_Depth"]
        result["n_ref_count"] = fmt_dict[normal_id][
            "Fragment_Supporting_Each_Alleles"
        ].split(",")[0]
        result["n_alt_count"] = fmt_dict[normal_id][
            "Fragment_Supporting_Each_Alleles"
        ].split(",")[-1]

    return result


def parse_sample_names(colnames: str) -> list[str]:
    """
    Recover sample names
    """
    return colnames.split("\t")[9:]


tumod_id: str | None = None
normal_id: str | None = None
parsed_lines: list[str] = []

with bz2.open(snakemake.input[0]) as bz2_stream:
    for line in bz2_stream:
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            sample_names = parse_sample_names(line[:-1])
            tumod_id = sample_names[0]
            if len(sample_names) > 1:
                normal_id = sample_names[-1]
        else:
            parsed_lines.append(
                line_to_maf(
                    line[:-1],
                    tumod_id=tumod_id,
                    headers_ANN=headers_ANN,
                    headers_fmt=headers_fmt,
                    normal_id=normal_id,
                )
            )

maf_table = pandas.DataFrame.from_records(parsed_lines)
maf_table.to_csv(
    snakemake.output[0],
    sep="\t",
    header=True,
    index=False,
)

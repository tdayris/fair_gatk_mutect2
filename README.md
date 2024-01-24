[![Snakemake](https://img.shields.io/badge/snakemake-≥7.29.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_gatk_mutect_germline/workflows/Tests/badge.svg?branch=main)](https://github.com/tdayris/fair_gatk_mutect_germline/actions?query=branch%3Amain+workflow%3ATests)

Do not use. Active dev.

Snakemake workflow used to call peaks with Macs2 and annotate them with Homer.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_gatk_mutect_germline) 
it is also available [locally](https://github.com/tdayris/fair_gatk_mutect_germline/blob/main/workflow/report/usage.rst) on a single page.
 
## Results

A complete description of the results can be found here in [workflow reports](https://github.com/tdayris/fair_gatk_mutect_germline/blob/main/workflow/report/results.rst).

## Material and Methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_gatk_mutect_germline/blob/main/workflow/report/material_methods.rst) textually. Web-links are available below:

![workflow_rulegraph](dag.png)

### Index and genome sequences with [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/tree/main)

#### Get DNA sequences

| Step                             | Commands                                                                                                         |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| Download DNA Fasta from Ensembl  | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-sequence.html) |
| Remove non-canonical chromosomes | [pyfaidx](https://github.com/mdshw5/pyfaidx)                                                                     |
| Index DNA sequence               | [samtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/faidx.html)                     |
| Creatse sequence Dictionary      | [picard](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/picard/createsequencedictionary.html)      |

#### Get genome annotation (GTF)

| Step                                                       | Commands                                                                                                             |
| ---------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| Download GTF annotation                                    | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-annotation.html) |
| Fix format errors                                          | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gff2gtf.html)                                     |
| Remove non-canonical chromosomes, based on above DNA Fasta | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sq_filter_feature_from_fasta.html)                           |
| Remove `<NA>` Transcript support levels                    | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sp_filter_feature_by_attribute_value.html)                   |

#### Get dbSNP variants

| Step                             | Commands                                                                                                                                     |
| -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| Download dbSNP variants          | [ensembl-variation](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-variation.html)                           |
| Filter non-canonical chromosomes | [pyfaidx](https://github.com/mdshw5/pyfaidx) + [BCFTools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bcftools/filter.html) |
| Index variants                   | [tabix](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/tabix/index.html)                                                       |

### Bowtie2 Mapping with [`fair_bowtie2_mapping`](https://github.com/tdayris/fair_bowtie2_mapping/tree/main)

#### Align reads over the genome

| Step             | Meta-Wrapper                                                                                                             | Wrapper                                                                                                                          |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| Bowtie2-build    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bowtie2/build.html)                                 |
| Fastp            |                                                                                                                          | [fastp](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastp.html)                                                 |
| Bowtie2-align    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bowtie2/align.html)                                 |
| Sambamba sort    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/sort.html)                                 |
| Sambamba-view    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/view.html)                                 |
| Sambamba-markdup | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/markdup.html)                           |
| Sambamba-index   | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-index](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/index.html)                               |

#### Quality controls

| Step     | Wrapper                                                                                                                          |
| -------- | -------------------------------------------------------------------------------------------------------------------------------- |
| Picard   | [picard-collectmultiplemetrics](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html) |
| Samtools | [samtools-stats](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/stats.html)                               |
| FastQC   | [fastqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/fastqc.html)                                       |
| MultiQC  | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/multiqc.html)                                     |


### Call variants with Mutect2

#### Actual Calling

| Step                                | Meta-Wrapper                                                                                                               | Wrapper                                                                                                                          |
| ----------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| Per-sample annotation               | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/gatk_mutect2_calling.html) | [add-or-replace-groups](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/picard/addorreplacereadgroups.html)         |
| Mutect2 calling                     | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/gatk_mutect2_calling.html) | [mutect2](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gatk/mutect.html)                                         |
| Infer contaminations                | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/gatk_mutect2_calling.html) | [get-pileup-summaries](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gatk/getpileupsummaries.html)                |
| Estimate corss-sample contamination | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/gatk_mutect2_calling.html) | [calculate-contamination](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gatk/calculatecontamination.html)         |
| Search for sequencing artifact bias | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/gatk_mutect2_calling.html) | [learn-read-orientation-model](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gatk/learnreadorientationmodel.html) |
| Filtering calls                     | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/gatk_mutect2_calling.html) | [filter-mutect-calls](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gatk/filtermutectcalls.html)                  |


#### Quality controls

| Step               | Wrapper                                                                                            |
| ------------------ | -------------------------------------------------------------------------------------------------- |
| Variant Evaluation | [variant-eval](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gatk/varianteval.html) |
| MultiQC            | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/multiqc.html)       |
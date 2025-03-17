[![Snakemake](https://img.shields.io/badge/snakemake-≥7.29.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_gatk_mutext2/workflows/Tests/badge.svg)](https://github.com/tdayris/fair_gatk_mutext2/actions?query=branch%3Amain+workflow%3ATests)


Snakemake workflow used to call germline variants with GATK-Mutect2

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_gatk_mutext2) 
it is also available [locally](https://github.com/tdayris/fair_gatk_mutext2/blob/main/workflow/report/usage.rst) on a single page.
 
## Results

A complete description of the results can be found here in [workflow reports](https://github.com/tdayris/fair_gatk_mutext2/blob/main/workflow/report/results.rst).

## Material and Methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_gatk_mutext2/blob/main/workflow/report/material_methods.rst) textually. Web-links are available below:

![workflow_rulegraph](dag.png)

### Index and genome sequences with [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/tree/main)

See [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/tree/main) for information about sequences and annotation retrieval.

### Raw-sequences QC with [`fair_fastqc_multiqc`](https://github.com/tdayris/fair_fastqc_multiqc/)

See  [`fair_fastqc_multiqc`](https://github.com/tdayris/fair_fastqc_multiqc/) documentation about ranw sequences quality controls

### Bowtie2 Mapping with [`fair_bowtie2_mapping`](https://github.com/tdayris/fair_bowtie2_mapping/tree/main)

See [`fair_bowtie2_mapping`](https://github.com/tdayris/fair_bowtie2_mapping/tree/main) for informatin about sequence alignment and quality controls.

### Call variants with Mutect2

#### Actual Calling

| Step                                | Meta-Wrapper                                                                                                               | Wrapper                                                                                                                          |
| ----------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| Per-sample annotation               | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/gatk_mutect2_calling.html) | [add-or-replace-groups](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/picard/addorreplacereadgroups.html)         |
| Mutect2 calling                     | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/gatk_mutect2_calling.html) | [mutect2](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/gatk/mutect.html)                                         |
| Infer contaminations                | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/gatk_mutect2_calling.html) | [get-pileup-summaries](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/gatk/getpileupsummaries.html)                |
| Estimate corss-sample contamination | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/gatk_mutect2_calling.html) | [calculate-contamination](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/gatk/calculatecontamination.html)         |
| Search for sequencing artifact bias | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/gatk_mutect2_calling.html) | [learn-read-orientation-model](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/gatk/learnreadorientationmodel.html) |
| Filtering calls                     | [GATK short variant calling](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/meta-wrappers/gatk_mutect2_calling.html) | [filter-mutect-calls](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/gatk/filtermutectcalls.html)                  |


```
┌─────────────────────────┐           ┌──────────────────────────────┐   
│Annotate samples (Picard)│           │Index annotated bam (Sambamba)│   
└─────────────┬───────────┘           └──────┬───────────────────────┘   
              │                              │                           
              ├──────────────────────────────┘                           
              │                                                          
┌─────────────▼─────────┐                                                
│Call variants (Mutect2)├──────────────────────┐                         
└─────────────┬─────────┘                      │                         
              │                                │                         
              │                                │                         
              │                                │                         
┌─────────────▼────────────┐         ┌─────────▼────────────────────────┐
│Infer contamination (GATK)│         │Estimate sequencing bias artifacts│
└─────────────┬────────────┘         │             (GATK)               │
              │                      └────────────┬─────────────────────┘
              │                                   │                      
┌─────────────▼─────────────────────┐             │                      
│Estimate cross-sample contamination│             │                      
│          (GATK)                   │             │                      
└─────────────────────────────────┬─┘             │                      
                                  │               │                      
                                  ├───────────────┘                      
                                  │                                      
                      ┌───────────▼──────────────────┐                   
                      │Hard filtering variants (GATK)│                   
                      └──────────────────────────────┘                   
```

#### Annotate variants of confidence

```
┌──────────────────────────┐     
│Annotate variants (SnpEff)│     
└──────────────┬───────────┘     
               │                 
               │                 
               │                 
┌──────────────▼────────────────┐
│Annotate variant type (SnpSift)│
└───────────────────────────────┘
```


#### Quality controls

| Step               | Wrapper                                                                                            |
| ------------------ | -------------------------------------------------------------------------------------------------- |
| Variant Evaluation | [variant-eval](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/gatk/varianteval.html) |
| MultiQC            | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v5.8.3/wrappers/multiqc.html)       |

```
┌───────────────────────────┐ ┌───────────────────┐ ┌────────────────────┐
│Variant annotation (SnpEff)│ │fair_fastqc_multiqc│ │fair_bowtie2_mapping│
└─────────────────────┬─────┘ └───┬───────────────┘ └───┬────────────────┘
                      │           │                     │                 
                      │           │                     │                 
                      └───────────┼─────────────────────┘                 
                                  │                                       
                                  │                                       
                     ┌────────────▼───────────────┐                       
                     │Report aggregation (MultiQC)│                       
                     └────────────────────────────┘                       
```

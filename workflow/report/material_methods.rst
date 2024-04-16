Matierial and methods
=====================

Genome DNA sequence and annotations were download from Ensembl. 
Pyfaidx_ [#pyfaidxpaper]_ was used to filter non-cannonical 
chromosomes. Agat_ [#agatpaper]_ was used to correct common 
issues found in Ensembl genome annotation files, filter non-
cannonical chromosomes, and remove transcripts with TSL being
equal to NA. Samtools_ [#samtoolspaper]_ and Picard_ [#picardpaper]_ 
were used to index genome sequences.

Raw fastq file quality was assessed with FastQC_ [#fastqcpaper]_.
Raw fastq files were trimmed using Fastp_ [#fastppaper]_ . Cleaned 
reads were aligned over indexed Ensembl genome with Bowtie2_ 
[#bowtie2paper]_. Sambamba_ [#sambambapaper]_ was used to sort, 
filter, mark duplicates, and compress aligned reads. Quality 
controls were done on cleaned, sorted, deduplicated aligned reads 
using Picard_ [#picardpaper]_ and Samtools_ [#samtoolspaper]_.
Goleft_ [#goleftpaper]_, RSeQC_ [#rseqcpaper]_ and ngsderive_ 
[#ngsderivepaper]_ were used to assess quality of mapping.


Read groups were re-defined over aligned reads with Picard_. GATK_
Mutect2_ [#mutect2paper]_ was used to perform germline calling,
following the `GATK Best practices`_ [#gatkbestpaper]_ and the
method described in the GATK germline best practices [#germlinepaper]_.
VariantEval_ [#gatkbestpaper]_ and BCFTools_ [#bcftoolspaper]_
were used to assess variant calling quality. All quality repors 
produced during both trimming and mapping steps have been aggregated 
with MultiQC_ [#multiqcpaper]_. 

The whole pipeline was powered by Snakemake_ [#snakemakepaper]_. 
This pipeline is freely available on Github_, details about 
installation usage, and resutls can be found on the 
`Snakemake workflow`_ page.

.. [#pyfaidxpaper] Shirley, Matthew D., et al. Efficient" pythonic" access to FASTA files using pyfaidx. No. e1196. PeerJ PrePrints, 2015.
.. [#agatpaper] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [#samtoolspaper] Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079.
.. [#picardpaper] McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303.
.. [#fastqcpaper] Andrews, S. Fastqc. "A quality control tool for high throughput sequence data. Augen, J.(2004). Bioinformatics in the post-genomic era: Genome, transcriptome, proteome, and information-based medicine." (2010).
.. [#fastppaper] Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics 34.17 (2018): i884-i890.
.. [#bowtie2paper] Langmead, Ben, and Steven L. Salzberg. "Fast gapped-read alignment with Bowtie 2." Nature methods 9.4 (2012): 357-359.
.. [#sambambapaper] Tarasov, Artem, et al. "Sambamba: fast processing of NGS alignment formats." Bioinformatics 31.12 (2015): 2032-2034.
.. [#goleftpaper] Pedersen, Brent S., et al. "Indexcov: fast coverage quality control for whole-genome sequencing." Gigascience 6.11 (2017): gix090.
.. [#rseqcpaper] Wang, Liguo, Shengqin Wang, and Wei Li. "RSeQC: quality control of RNA-seq experiments." Bioinformatics 28.16 (2012): 2184-2185.
.. [#ngsderive] McLeod, Clay, et al. "St. Jude Cloud: a pediatric cancer genomic data-sharing ecosystem." Cancer discovery 11.5 (2021): 1082-1099.
.. [#mutect2paper] Benjamin, David, et al. "Calling somatic SNVs and indels with Mutect2." BioRxiv (2019): 861054.
.. [#gatkbestpaper] Van der Auwera, Geraldine A., and Brian D. O'Connor. Genomics in the cloud: using Docker, GATK, and WDL in Terra. O'Reilly Media, 2020.
.. [#germlinepaper] Poplin, Ryan, et al. "Scaling accurate genetic variant discovery to tens of thousands of samples." BioRxiv (2017): 201178.
.. [#bcftoolspaper] Danecek, Petr, et al. "Twelve years of SAMtools and BCFtools." Gigascience 10.2 (2021): giab008.
.. [#multiqcpaper] Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.

.. _Sambamba: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/sambamba.html
.. _Bowtie2: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/bowtie2.html
.. _Fastp: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/fastp.html
.. _Picard: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/picard/collectmultiplemetrics.html
.. _MultiQC: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/multiqc.html
.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_bowtie2_mapping
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_bowtie2_mapping
.. _Agat: https://agat.readthedocs.io/en/latest/index.html
.. _Samtools: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/samtools/faidx.html
.. _FastQC: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/fastqc.html
.. _Pyfaidx: https://github.com/mdshw5/pyfaidx
.. _GATK: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/gatk.html
.. _`GATK Best practices`: https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices
.. _Mutect2: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/gatk/mutect.html
.. _VariantEval: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/gatk/varianteval.html
.. _BCFTools: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/bcftools/stats.html
.. _ngsderive: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/sambamba.html
.. _rseqc: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/sambamba.html
.. _goleft: https://snakemake-wrappers.readthedocs.io/en/v3.7.0/wrappers/sambamba.html

:Authors:
    Thibault Dayris

:Version: 1.4.0 of 04/03/2024

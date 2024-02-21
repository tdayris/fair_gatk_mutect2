Results
=======


Alongside with the report, you may find directories called `reference`,
and `results`.


Results list
============


Reference
---------


Alongside with this report, you may find a directory called `reference`.
You shall find all requested files in it. By default, the following
files are present:

::

    reference/
    ├── blacklist
    |   └── XXX.merged.bed
    ├── variants
    |   ├── XXX.all.vcf.gz
    |   └── XXX.all.vcf.gz.tbi
    ├── sequences
    |   ├── XXX.cdna.fasta
    |   ├── XXX.cdna.fasta.fai
    |   ├── XXX.dna.dict
    |   ├── XXX.dna.fasta
    |   └── XXX.dna.fasta.fai
    └── annotation
        ├── XXX.id_to_gene.tsv
        ├── XXX.t2g.tsv
        └── XXX.gtf


+-------------------+-----------------------------+
| Extension         | Content                     |
+===================+=============================+
| `.bed`            | Genome blacklisted regions  |
+-------------------+-----------------------------+
| `.gtf`            | Genome annotation           |
+-------------------+-----------------------------+
| `.id_to_gene.tsv` | Genome id-to-name           |
+-------------------+-----------------------------+
| `.t2g.tsv`        | Transcript id-to-name       |
+-------------------+-----------------------------+
| `.fasta`          | Genome sequences            |
+-------------------+-----------------------------+
| `.fasta.fai`      | Genome sequences index      |
+-------------------+-----------------------------+
| `.dict`           | Genome sequences dictionary |
+-------------------+-----------------------------+
| `.vcf.gz`         | Genome known variations     |
+-------------------+-----------------------------+

These files are quite volumous and are not embeded in this HTML page. Please
find them directly on file system.


Results
-------

Given a samples called `YYY` and a genome called `XXX`,
the following files are present:


::

    results/
    ├── QC
    │   ├── MultiQC_FastQC_data.zip
    │   ├── MultiQC_FastQC.html
    │   ├── report_pe
    │   │   ├── YYY.1_fastqc.zip
    │   │   ├── YYY.1.html
    │   │   ├── YYY.2_fastqc.zip
    │   │   ├── YYY.2.html
    │   │   └── YYY.html
    └── XXX.dna
        ├── Mapping
        │   ├── YYY.bam
        │   ├── YYY.bam.bai
        ├── QC
        │   ├── MultiQC_GatkGermlineCalling_data.zip
        │   ├── MultiQC_GatkGermlineCalling.html
        │   ├── MultiQC_Mapping_data.zip
        │   └── MultiQC_Mapping.html
        └── VariantCalling
            └── Germline
                ├── YYY.vcf.gz
                ├── YYY.vcf.gz.filteringStats.tsv
                └── YYY.vcf.gz.tbi




+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
| Directory                | File Extension                          | Content                                                                      |
+==========================+=========================================+==============================================================================+
| XXX/Mapping              | `YYY.bam`                               | Aligned reads                                                                |
+                          +-----------------------------------------+------------------------------------------------------------------------------+
|                          | `YYY.bam.bai`                           | Aligned reads index                                                          |
+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
| XXX/VariantCalling/Raw   | `YYY.vcf.gz`                            | Germline variants, unannotated, unfiltered                                   |
+                          +-----------------------------------------+------------------------------------------------------------------------------+
|                          | `YYY.vcf.gz.tbi`                        | Germline variants index                                                      |
+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
| QC                       | `MultiQC_FastQC_data.zip`               | Zipped figures and tables                                                    |
+                          +-----------------------------------------+------------------------------------------------------------------------------+
|                          | `MultiQC_FastQC.html`                   | Complete quality report, includes all samples                                |
+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
|                          | `MultiQC_Mapping_data.zip`              | Zipped figures and tables                                                    |
+                          +-----------------------------------------+------------------------------------------------------------------------------+
|                          | `MultiQC_Mapping.html`                  | Complete quality report, includes all samples up to mapping step             |
+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
|                          | `MultiQC_GatkGermlineCalling_data.zip`  | Zipped figures and tables                                                    |
+                          +-----------------------------------------+------------------------------------------------------------------------------+
|                          | `MultiQC_GatkGermlineCalling.html`      | Complete quality report, includes all samples up to germline variant calling |
+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
| QC/report_pe             | `YYY.html`                              | Sequence quality report for PE sample `YYY`                                  |
+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
| QC/report_se             | `YYY.html`                              | Sequence quality report for SE sample `YYY`                                  |
+--------------------------+-----------------------------------------+------------------------------------------------------------------------------+
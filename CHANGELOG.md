# 1.5.5

## Features:

* Possible alignment with STAR
* fair_bowtie2_mapping up to 4.4.6
* fair_fastqc_multiqc up to 2.2.5
* fair_genome_indexer up to 3.9.7
* fair_star_mapping up to 1.3.3


# 1.5.4

## Features

* fair_bowtie2_mapping up to 4.1.3

# 1.5.3

## Features:

* fair_bowtie2_mapping up to 4.1.2

# 1.5.2

## Features:

* Better resources reservation
* snakemake-wrappers up to 3.13.7
* fair_bowtie2_mapping up to 4.1.0
* Allow local pipeline search
* Allow local snakemake-wrappers

## Fixes:

* Format error

# 1.5.1

## Features:

* fair_bowtie2_mapping up to 3.5.2
* fair_genome_indexer up to 3.8.1
* Memory / time reservations

# 1.5.0

## Features:

* Split multi-allelic sites and left-normalize
* Automatically filter passing variants
* Support of tumor/normal-pairs calling

# 1.4.2

## Features:

* More time for Mutect2 rule

# 1.4.1

## Features:

* fair_bowtie2_mapping update to 3.3.3
* fair_fastqc_multiqc update to 2.2.7

# 1.4.0

## Features:

* Snakemake wrappers up to 3.7.0
* fair_fastqc_multiqc up to version 2.2.6
* fair_bowtie2_mapping up to version 3.3.2
* fair_genome_indexer up to version 3.4.4


# 1.3.0

## Features:

* Handle somatic calling if possible
* Snakemake wrappers up to version 3.5.2
* fair_bowtie2_mapping up to version 3.3.0
* fair_fastqc_multiqc up to version 2.2.1
* fair_genome_indexer up to version 3.3.0
* Configuration keys are *all* optional

# 1.2.0

## Features:

* Snakemake wrappers up to version 3.5.0
* fair_bowtie2_mapping up to version 3.2.0
* fair_fastqc_multiqc up to version 2.1.2
* fair_genome_indexer up to version 3.2.2
* SnpEff + SnpSift basic annotations
* Pipeline containerized

# 1.1.1

## Fix:

* missing Log/benchmark

# 1.1.0

## Features:

* temp, log, and benchmark paths rebuild
* rule names changed to follow between-workflows trace
* use of `lookup` instead of hand-made function
* Relies on fair_genome_indexer version 3.1.4
* Relies on fair_bowtie2_mapping version 3.1.0
* Relies on fair_fastqc_multiqc version 2.0.3
* MultiQC now holds only one `{species}.{build}.{release}`

## Fixes:

* Documentation update

# 1.0.0

## Features:

* Snakemake v8+ compatible
* Snakemake-wrappers
* Relies on fair_genome_indexer version 3.0.0
* Relies on fair_bowtie2_mapping version 3.0.1
* Deployable workflow with integration tests on Github

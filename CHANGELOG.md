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

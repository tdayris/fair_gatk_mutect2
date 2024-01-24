This pipeline requires two configuration file:

# `config.yaml`

A standard `Snakemake` configuration, yaml-formatted file containing a list of
all parameters accepted in this workflow:

* `samples`: Path to the file containing link between samples and their fastq file(s)
* `params`: Per-tool list of optional parameters

Example:

```
samples: config/samples.csv

# Optional parameters
params:
  # Optional parameters for pyfaidx (filter/correct fasta format)
  pyfaidx:
    # Filter-out non canonical chromosomes
    dna: '--regex "^[[0-9]+|X|Y|MT]"'
    # Keep all cdna sequences
    cdna: ""
  # Optional parameters for agat (filter/correct GTF format)
  agat:
    # Optional parameters for agat_convert_sp_gff2gtf.pl
    gff2gtf: ""
    # Optional parameters for agat_sq_filter_feature_from_fasta.pl
    filter_features: ""
    # Optional parameters for agat_sq_select_feature_by_attribute_value.pl
    select_feature_by_attribute_value: "--attribute 'transcript_support_level' --value '\"NA\"' --test '='"
    # Optional parameters for agat_convert_sp_gff2tsv
    agat_convert_sp_gff2tsv: ""
  # Optional parameters for GFFRead
  gffread: ""
  # Optional parameters for bedtools
  bedtools:
    # Optional parameters for filtering non-canonical chromosomes over dbSNP
    filter_non_canonical_chrom: ""
  # Optional parameters for tabix index VCF
  tabix: "-p vcf"
  # Optional parameters for fastp
  fastp:
    # Optional adapters to remove
    adapters: ""
    # Optional command line arguments for fastp
    extra: ""
  # Optional parameters for fastqc
  fastqc: ""
  # Optional parameters for bowtie2
  bowtie2:
    # Optional parameters for bowtie2-build
    build: ""
    # Optional parameters for bowtie2-align
    align: ""
  sambamba:
    # Optional parameters for sambamba view
    view: "--format 'bam' --filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)' "
    # Optional parameters for sambamba markdup
    markdup: "--remove-duplicates --overflow-list-size=500000"
  picard:
    # Mapping QC optional parameters
    metrics: ""
    # Optional parameters for picard create sequence dictionary
    createsequencedictionary: ""
  # Optional parameters for samtools stats
  samtools:
    # Optional parameters for samtools fasta index
    faidx: ""
    # Optional parameters for samtools stats
    stats: ""
  # Optional parameters for multiqc
  multiqc: "--module picard --module fastqc --module fastp --module samtools --module bowtie2 --module sambamba --zip-data-dir --verbose --no-megaqc-upload --no-ansi --force"
  # Optional parameters for GATK
  gatk:
    # Optional parameters for Mutect2
    mutect2: ""
    # Optional parameters for GATK getpileupsummaries
    getpileupsummaries: ""
    # Optional parameters for GATK calculatecontamination
    calculatecontamination: ""
    # Optional parameters for GATK learnreadorientationmodel
    learnreadorientationmodel: ""
    # Optional parameters for GATK filtermutectcalls
    filtermutectcalls: "--create-output-variant-index --min-median-mapping-quality 35 --max-alt-allele-count 3"
    # Optional parameters for GATK varianteval
    varianteval: ""
```

# `samples.csv`

A CSV-formatted text file containing the following mandatory columns:

* sample_id: Unique name of the sample
* upstream_file: Path to upstream fastq file
* species: The species name, according to Ensembl standards
* build: The corresponding genome build, according to Ensembl standards
* release: The corresponding genome release, according to Ensembl standards
* downstream_file: Optional path to downstream fastq file

Example:

```
sample_id,upstream_file,downstream_file,species,build,release
sac_a,data/reads/a.scerevisiae.1.fq,data/reads/a.scerevisiae.2.fq,saccharomyces_cerevisiae,R64-1-1,110
sac_a_input,data/reads/a.scerevisiaeI.1.fq,data/reads/a.scerevisiaeI.2.fq,saccharomyces_cerevisiae,R64-1-1,110
```

While `CSV` format is tested and recommended, this workflow uses python
`csv.Sniffer()` to detect column separator. Tabulation and semicolumn are
also accepted as field separator. Remember that only comma-separator is
tested.
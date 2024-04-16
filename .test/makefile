SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt report.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/bcftools.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/rules/fair_bowtie2_mapping_pipeline.smk" \
	"../workflow/rules/fair_fastqc_multiqc_pipeline.smk" \
	"../workflow/rules/fair_genome_indexer_pipeline.smk" \
	"../workflow/rules/gatk_mutect2_calling_meta.smk" \
	"../workflow/rules/gatk.smk" \
	"../workflow/rules/multiqc.smk" \
	"../workflow/rules/snpeff.smk" \
	"../workflow/rules/snpsift.smk" \
	"../workflow/Snakefile"


# List of python files to lint
python_scripts := \
	"../workflow/scripts/fair_gatk_mutect_germline_multiqc_config.py" 



.PHONY: all
all: report.txt


dag.png: format.txt linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--rulegraph --forceall | dot -Tpng > dag.png && \
	cp --verbose dag.png ../dag.png


report.txt: pipeline.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --profile local \
		--report report.zip 2>&1 | tee report.txt


pipeline.txt: linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--cores 7 \
		--software-deployment-method conda \
		2>&1 | tee pipeline.txt


format.txt:
	mamba activate test-snakemake-wrappers && \
	snakefmt $(snakefiles) 2>&1 | tee format.txt


linter.txt: format.txt black.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--lint 2>&1 | tee linter.txt


black.txt:
	mamba activate test-snakemake-wrappers && \
	black $(python_scripts) 2>&1 | tee black.txt


.PHONY: clean
clean:
	rm -rvf linter.txt \
		   format.txt \
		   reference \
		   results \
		   tmp \
		   logs \
		   benchmark \
		   report.html \
		   report.txt \
		   pipeline.txt \
		   black.txt \
		   report.zip \
		   report

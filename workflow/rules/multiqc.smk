rule fair_gatk_mutect_germline_multiqc_report:
    input:
        unpack(get_fair_gatk_mutect_germline_multiqc_report_input),
    output:
        report(
            "results/QC/MultiQC_GatkGermlineCalling.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "GatkGermlineCalling",
            },
        ),
        "results/QC/MultiQC_GatkGermlineCalling_data.zip",
    params:
        extra=config.get("params", {}).get(
            "multiqc", "--zip-data-dir --verbose --no-megaqc-upload --no-ansi --force"
        ),
        use_input_files_only=True,
    log:
        "logs/fair_gatk_mutect_germline/multiqc_report.log",
    benchmark:
        "benchmark/fair_gatk_mutect_germline/multiqc_report.tsv"
    wrapper:
        "v3.3.6/bio/multiqc"

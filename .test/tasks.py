# coding: utf-8

from invoke import task
from github import Github
from datetime import datetime

import os.path as op
import os
import yaml

tag_re = r"v\?[0-9]\+\.[0-9]\+\.[0-9]\+"

locations = {
    "snakefile": "../workflow/Snakefile",
    "rules": "../workflow/rules/",
    "scripts": "../workflow/scripts",
    "envs": "../workflow/envs",
    "changelog": "../CHANGELOG.md",
    "cff": "../CITATION.cff",
}
for name, location in locations.items():
    if not op.exists(location):
        raise FileNotFoundError(f"Could not find {location=}")

    else:
        locations[name] = op.realpath(location)


def get_future_version(
    changelog: str = locations["changelog"],
) -> str:
    with open(changelog, "r") as changelog_stream:
        line = next(changelog_stream)
        version = line.strip("#").strip()

    print(f"Future version is: {version=}")
    return version


future_version = get_future_version()


def get_latest_release(
    address: str,
    future: str = future_version,
    known: str | None = None,
) -> str:
    name = address.split("/")[-1]
    current_name = op.basename(op.dirname(os.getcwd()))
    if name == current_name:
        print(f"Current pipeline detected, using {future=} for {name=}...")
        return get_future_version(changelog=locations["changelog"])

    if known is not None:
        print(f"Known version provided {known=} for {name=}")
        return known

    print(f"Seaching for {address=} on github...")
    git = Github()
    repo = git.get_repo(address)
    releases = repo.get_releases()
    return releases[0].tag_name

snakemake_wrappers_version = get_latest_release(
    address="snakemake/snakemake-wrappers",
    future=future_version,
    known=os.environ.get("SNAKEMAKE_WRAPPERS_VERSION"),
)

fair_genome_indexer_version = get_latest_release(
    "tdayris/fair_genome_indexer",
    future=future_version,
    known=os.environ.get("FAIR_GENOME_INDEXER_VERSION"),
)

fair_fastqc_multiqc_version = get_latest_release(
    "tdayris/fair_fastqc_multiqc",
    future=future_version,
    known=os.environ.get("FAIR_FASTQC_MULTIQC_VERSION"),
)

fair_bowtie2_mapping_version = get_latest_release(
    "tdayris/fair_bowtie2_mapping",
    future=future_version,
    known=os.environ.get("FAIR_BOWTIE2_MAPPING_VERSION"),
)

fair_star_mapping_version = get_latest_release(
    "tdayris/fair_star_mapping",
    future=future_version,
    known=os.environ.get("FAIR_STAR_MAPPING_VERSION"),
)

fair_macs2_calling_version = get_latest_release(
    "tdayris/fair_macs2_calling",
    future=future_version,
    known=os.environ.get("FAIR_MACS2_CALLING_VERSION"),
)

fair_gatk_mutect2_calling_version = get_latest_release(
    "tdayris/fair_gatk_mutect2_calling",
    future=future_version,
    known=os.environ.get("FAIR_GATK_MUTECT2_CALLING_VERSION"),
)

@task
def clean(
    c,
    conda: bool = False,
    extra: str = "",
):
    patterns = [
        "black.txt",
        "format.txt",
        "linter_info.txt",
        "pipeline.txt",
        "results",
        "resources",
        "report.txt",
        "report.zip",
        "resources.md",
        "resources.tsv",
        "resources.txt",
        "summary.tsv",
        "summary_small.tsv",
        "docs_update.txt",
        "logs",
        "reference",
        "report",
        "tmp",
        "Snakefile",
        "wrappers_update.txt",
        "conda_update.txt",
        "update_docs.txt",
        "docs_update.txt",
        "linter.txt",
    ]
    if conda:
        patterns.append(".snakemake")
        patterns.append(".conda")

    for pattern in patterns:
        if op.exists(pattern):
            c.run(
                f"rm --force --recursive --verbose '{pattern}'",
                echo=True,
            )


@task
def update_docs_cff(
    c,
    to: str = future_version,
    cff_path: str = locations["cff"],
):
    today: str = datetime.today().strftime("%Y-%m-%d")
    name: str = op.basename(op.dirname(os.getcwd()))
    cff = {
        "cff-version": "1.2.0",
        "message": "If you use this software, please cite it as below.",
        "authors": [
            {
                "family-names": "Dayris",
                "given-names": "Thibault",
                "orcid": "https://orcid.org/0009-0009-2758-8450",
            }
        ],
        "title": name.replace("_", "-"),
        "version": future_version,
        "date-released": today,
        "url": f"https://github.com/tdayris/{name}",
    }
    print(cff)
    with open(cff_path, "w") as yaml_cff_stream:
        yaml.dump(cff, yaml_cff_stream, default_flow_style=False)


@task(update_docs_cff)
def update_docs_wrappers(
    c,
    to: str = snakemake_wrappers_version,
    future_version: str = future_version,
    tag_re: str = tag_re,
):
    today = datetime.today().strftime("%Y-%m-%d")
    regex = rf"s|{tag_re}/wrappers|{to}/wrappers|g"
    print(f"Updating snakemake wrappers in README.md to {to=}")
    c.run(
        f"sed -i '{regex}' ../README.md >> docs_update.txt 2>&1",
        echo=True,
    )

    for root, dirs, files in os.walk("../workflow/report"):
        for file in files:
            if file.endswith(".rst"):
                print(f"Updating snakemake wrappers in '{root}/{file}'...")
                c.run(
                    f"sed -i '{regex}' '{root}/{file}' >> docs_update.txt 2>&1",
                    echo=True,
                )

    regex = (
        's|snakemake_wrappers_prefix: str = "v4.5.0"|'
        f'snakemake_wrappers_prefix: str = "{to}"|g'
    )
    print("Updating '../workflow/rules/common.smk'...")
    c.run(
        f"sed -i '{regex}' '../workflow/rules/common.smk' >> update_docs.txt 2>&1",
        echo=True,
    )

    regex = (
        rf"s|fair-genome-indexer (Version {tag_re})"
        f"|fair-genome-indexer (Version {fair_genome_indexer_version})|g;"
        rf"s|fair-fastqc-multiqc (Version {tag_re})"
        f"|fair-fastqc-multiqc (Version {fair_fastqc_multiqc_version})|g;"
        rf"s|fair-bowtie2-mapping (Version {tag_re})"
        f"|fair-bowtie2-mapping (Version {fair_bowtie2_mapping_version})|g;"
        rf"s|fair-star-mapping (Version {tag_re})"
        f"|fair-star-mapping (Version {fair_star_mapping_version})|g;"
        rf"s|fair-macs2-calling (Version {tag_re})"
        f"|fair-macs2-calling (Verison {fair_macs2_calling_version})|g;"
        rf"s|fair-gark-mutect2 (Version {tag_re})"
        f"|fair-gatk-mutect2 (Verison {fair_gatk_mutect2_calling_version})|g"
    )
    print("Updating '../workflow/report/material_methods.rst'...")
    c.run(
        f"sed -i '{regex}' '../workflow/report/material_methods.rst' >> update_docs.txt 2>&1",
        echo=True,
    )

    regex = (
        rf"s|\:Version\: {tag_re} of [0-9]\+\-[0-9]\+\-[0-9]\+$"
        rf"|\:Version\: {future_version} of {today}|g"
    )
    c.run(
        f"sed -i '{regex}' '../workflow/report/material_methods.rst' "
        ">> update_docs.txt 2>&1",
        echo=True,
    )


@task(update_docs_wrappers)
def update_wrappers_rules(
    c,
    to: str = snakemake_wrappers_version,
    snakefile: str = locations["snakefile"],
    rules: str = locations["rules"],
    tag_re: str = tag_re,
):
    print(f"Updating {snakefile=}...")
    c.run(
        "snakedeploy update-snakemake-wrappers "
        f"--git-ref '{snakemake_wrappers_version}' '{snakefile}' "
        ">> wrappers_update.txt 2>&1",
        echo=True,
    )

    for root, dirs, files in os.walk(rules):
        for file in files:
            if file == "fair_genome_indexer.smk":
                c.run(
                    rf"""sed -i 's|tag="{tag_re}"|tag="{fair_genome_indexer_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""",
                    echo=True,
                )
            elif file == "fair_fastqc_multiqc.smk":
                c.run(
                    rf"""sed -i 's|tag="{tag_re}"|tag="{fair_fastqc_multiqc_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""",
                    echo=True,
                )
            elif file == "fair_bowtie2_mapping.smk":
                c.run(
                    rf"""sed -i 's|tag="{tag_re}"|tag="{fair_bowtie2_mapping_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""",
                    echo=True,
                )
            elif file == "fair_star_mapping.smk":
                c.run(
                    rf"""sed -i 's|tag="{tag_re}"|tag="{fair_star_mapping_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""",
                    echo=True,
                )
            elif file == "fair_macs2_calling.smk":
                c.run(
                    rf"""sed -i 's|tag="{tag_re}"|tag="{fair_macs2_calling_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""",
                    echo=True,
                )
            elif file == "fair_gatk_mutect2.smk":
                c.run(
                    rf"""sed -i 's|tag="{tag_re}"|tag="{fair_gatk_mutect2_calling_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""",
                    echo=True,
                )
            elif file.endswith(".smk"):
                c.run(
                    "snakedeploy update-snakemake-wrappers --git-ref "
                    f"'{snakemake_wrappers_version}' '{root}/{file}' "
                    ">> wrappers_update.txt 2>&1",
                    echo=True,
                )


@task(update_wrappers_rules)
def update_conda(
    c,
    envs: str = locations["envs"],
):
    for root, dirs, files in os.walk(envs):
        for file in files:
            if file.endswith(".yaml"):
                c.run(
                    "snakedeploy update-conda-envs --conda-frontend mamba "
                    f"--pin-envs '{root}/{file}' > conda_update.txt 2>&1",
                    echo=True,
                )


@task
def black(c, scripts: str = locations["scripts"]):
    log: str = "black.txt"
    for root, dirs, files in os.walk(scripts):
        for file in files:
            if file.endswith(".py"):
                c.run(
                    f"black '{root}/{file}' >> {log} 2>&1",
                    echo=True,
                )


@task(black)
def snakefmt(
    c,
    snakefile: str = locations["snakefile"],
    rules: str = locations["rules"],
):
    log: str = "format.txt"
    c.run(
        f"snakefmt '{snakefile}' >> '{log}' 2>&1",
        echo=True,
    )

    for root, dirs, files in os.walk(rules):
        for file in files:
            if file.endswith(".smk"):
                c.run(
                    f"snakefmt '{root}/{file}' >> '{log}' 2>&1",
                    echo=True,
                )


@task(snakefmt)
def linter(
    c,
    snakefile: str = locations["snakefile"],
):
    c.run(
        f"snakemake --lint -s '{snakefile}' > linter_info.txt 2>&1",
        echo=True,
    )


@task(linter)
def pipeline(
    c,
    snakefile: str = locations["snakefile"],
):
    c.run(
        f"snakemake -s '{snakefile}' "
        "--cores 7 --restart-times 0 "
        "--rerun-incomplete --printshellcmds "
        "--shadow-prefix 'tmp' --rerun-triggers 'mtime' "
        "--software-deployment-method conda "
        "--benchmark-extended > pipeline.txt 2>&1",
        echo=True,
    )


@task(pipeline)
def report(
    c,
    snakefile: str = locations["snakefile"],
):
    c.run(
        f"snakemake -s '{snakefile}' --report report.zip > report.txt 2>&1",
        echo=True,
    )


@task(clean, update_conda, report)
def all(c):
    print("All done.")

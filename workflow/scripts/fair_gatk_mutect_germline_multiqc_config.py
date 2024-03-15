# -*- coding: utf-8 -*-

"""Snakemake wrapper for MultiQC configuration file"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import yaml

with open(str(snakemake.output[0]), "w") as out_yaml_stream:
    out_yaml_stream.write(
        yaml.dump(snakemake.params["extra"], default_flow_style=False)
    )

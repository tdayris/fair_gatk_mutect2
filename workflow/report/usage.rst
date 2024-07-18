`Gustave Roussy`_ special launching steps
=========================================

If you belong to `Gustave Roussy`_, and use `flamingo` computing cluster, please use the following command lines
and ignore the rest of this documentation.

::

    # Activate conda environment
    # Use latest version available, eg:
    conda activate /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_install/snakemake

    # Deploy workflow with the version of your choice
    snakedeploy deploy-workflow \
        https://github.com/tdayris/fair_gatk_mutect_germline . \
        --tag <version>

    # Run snakemake command
    snakemake --profile '/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/profiles/slurm-web/'

With `version` being the latest available version of this pipeline. Select your `version here`_


Step 1 : Install Snakemake and Snakedeploy
==========================================

Snakemake_ and Snakedeploy_ are best installed via the `Mamba package manager`_ 
(a drop-in replacement for conda). If you have neither Conda nor Mamba, it can 
be installed via Mambaforge_. For other options see: `mamba-org`_.

Given that Mamba is installed, run

::
    
    mamba create -c conda-forge \
                 -c bioconda \
                 --name snakemake \
                 snakemake \
                 snakedeploy \
                 mamba \


to install both Snakemake and Snakedeploy in an isolated environment.
For all following commands ensure that this environment is activated 
via the following command:

::
    
    conda activate snakemake


Step 2 : Deploy workflow
========================

Given that Snakemake_ and Snakedeploy_ are installed and available (see Step 1),
the workflow can be deployed as follows.

First, create an appropriate project working directory on your system and enter it:

::
    
    mkdir -p path/to/project-workdir
    cd path/to/project-workdir

In all following steps, we will assume that you are inside of that directory.

Second, run:

::
    
    snakedeploy deploy-workflow \
                https://github.com/tdayris/fair_gatk_mutect_germline . \
                --tag <version>

Where <version> is the latest available verison.

Snakedeploy will create two folders `workflow` and `config`. The former contains the 
deployment of the chosen workflow as a `Snakemake module`_, the latter contains 
configuration files which will be modified in the next step in order to configure 
the workflow to your needs. Later, when executing the workflow, Snakemake will 
automatically find the main `Snakefile` in the `workflow` subfolder.

Third, consider to put this directory under version control, e.g. by 
`managing it via a (private) Github repository`_


Step 3 : Configure the workflow
===============================

Edit the file `config.yaml` and `genomes.csv` according to the description
available in the `config/README.md`_ file.

Step 4: Run workflow
====================

Given that the workflow has been properly deployed and configured, it can be executed 
as follows.

Fow running the workflow while deploying any necessary software via conda (using 
the `Mamba package manager`_ by default), run Snakemake with:

::
    
    snakemake --cores all --software-deployment-method conda

Snakemake will automatically detect the main `Snakefile` in the `workflow` subfolder 
and execute the workflow module that has been defined by the deployment in step 2.

For further options, e.g. for cluster and cloud execution, see Snakemake_ documentation.

Step 5 : Generate report
========================

After finalizing your data analysis, you can automatically generate an interactive visual 
HTML report for inspection of results together with parameters and code inside of the 
browser via:

::
    
    snakemake --report report.zip

The resulting `report.zip` file can be passed on to collaborators, provided as a supplementary 
file in publications, or uploaded to a service like Zenodo_ in order to obtain a citable DOI_. 

.. _Snakemake: https://snakemake.readthedocs.io/en/stable/index.html
.. _Snakedeploy: https://snakedeploy.readthedocs.io/en/latest/
.. _`Mamba package manager`: https://github.com/mamba-org/mamba
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _`mamba-org`: https://github.com/mamba-org/mamba
.. _`Snakemake module`: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows
.. _`managing it via a (private) Github repository`: https://docs.github.com/en/github/importing-your-projects-to-github/adding-an-existing-project-to-github-using-the-command-line
.. _`config/README.md`: https://github.com/tdayris/fair_gatk_mutect_germline/blob/main/config/README.md
.. _Zenodo: https://zenodo.org/
.. _DOI: https://en.wikipedia.org/wiki/Digital_object_identifier
.. _`Gustave Roussy`: https://www.gustaveroussy.fr/en
.. _`version here`: https://github.com/tdayris/fair_gatk_mutect_germline/releases

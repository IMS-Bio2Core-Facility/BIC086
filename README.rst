BIC086-Sophie-Austin
====================

GTEx isoform expression queries for transcription factors in human hypothalamus.

Motivation
----------

When choosing transcription factors (TFs) for over expression iPSCs,
we stand the highest chance of success by choosing open reading frames
that are expressed in our cell-type or tissue of interest.
For tissue-dependent expression data,
there are few resources better than GTEx.
In this case, the `medianTranscriptExpression` query provides the necessary data.
It returns the median expression of each transcript for a gene in a given tissue.
Here, we query a list of genes against the 202 brain hypothalamus datasets.

Data
----

The pipeline requires no input data other than a list of gene names specified input
`configuration/snakemake.yaml`.

Reference
---------

Under the hood, GTEx uses Gencode v26 for its Ensembl IDs.
As this is not the most up-to-date version,
it actually proved quite frustrating to find the desired version numbers for each gene.
To streamline this process,
the pipeline now downloads the Gencode v26 GTF annotations from EBI,
removes all features not annotated as `gene`,
and pulls the names and IDs of all remaining features.
This way, all end users need to do is provide a list of gene names.

Pipeline
--------

Installation
~~~~~~~~~~~~

This pipeline needs `conda`_ and `snakemake`_ installed,
and runs best if you also have `singularity`_ installed.

Snakemake recommends using `mambaforge`_ as your base conda,
which I would also recommend.
Installation instructions are at the above link.
If you prefer a vanilla conda installation,
you can always try `mamba` following the instructions at the above snakemake link.
Once you have conda installed,
install snakemake as outlined on their page
(again, see the above link).

Singularity can be a bit of a nuisance to install
(at least in my experience).
I would recommend a Linux OS
(it's most straightforward on CentOS/RHEL  - hello AWS!)
and following the instructions in their repository `here`_.

.. _conda: https://docs.conda.io/en/latest/
.. _snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
.. _singularity: https://sylabs.io/singularity/
.. _mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _here: https://github.com/sylabs/singularity/blob/master/INSTALL.md


Running
~~~~~~~

Once you've installed the above software,
running the pipeline is as simple as:

.. code:: shell

   git clone https://github.com/IMS-Bio2Core-Facility/BIC086-Sophie-Austin &&
   cd BIC086-Sophie-Austin &&
   snakemake --use-conda --use-singularity --cores 6

Unless you are querying a huge number of genes,
I find 6 cores to be sufficient to keep things moving.
If you want to use more or less,
alter the `threads` parameter in `configuration/snakemake.yaml`.
That file is also the same spot to specify your gene list of interest!

If you aren't using `singularity`,
then leave off the apropriate flag, as so:

.. code:: shell

   git clone https://github.com/IMS-Bio2Core-Facility/BIC086-Sophie-Austin &&
   cd BIC086-Sophie-Austin &&
   snakemake --use-conda --cores 6

Under

BIC086-Sophie-Austin
====================

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: MIT License

.. image:: https://img.shields.io/badge/Python-3.9-brightgreen.svg
   :target: https://docs.python.org/3/whatsnew/3.9.html
   :alt: Python 3.9

.. image:: https://www.repostatus.org/badges/latest/active.svg
   :alt: Project Status: Active – The project has reached a stable, usable state and is being actively developed.
   :target: https://www.repostatus.org/#active

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Codestyle: Black

.. image:: https://img.shields.io/github/stars/IMS-Bio2Core-Facility/BIC086-Sophie-Austin?style=social
   :target: https://github.com/IMS-Bio2Core-Facility/BIC086-Sophie-Austin
   :alt: GitHub Repo stars

A fully concurrent pipeline for querying transcript-level GTEx data in the human hypothalamus.

Motivation
----------

When choosing transcription factors (TFs) for over expression iPSCs,
we stand the highest chance of success by choosing open reading frames
that are expressed in our cell-type or tissue of interest.
For tissue-dependent expression data,
there are few resources better than GTEx.
In this case, the ``medianTranscriptExpression`` query provides the necessary data.
It returns the median expression of each transcript for a gene in a given tissue.
Here, we query a list of genes against the 202 brain hypothalamus datasets.

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
you can always try ``mamba`` following the instructions at the above snakemake link.
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

.. code-block:: shell

   git clone https://github.com/IMS-Bio2Core-Facility/BIC086-Sophie-Austin &&
   cd BIC086-Sophie-Austin &&
   snakemake --use-conda --use-singularity --cores 6

If you aren't using ``singularity``,
then leave off the apropriate flag, as so:

.. code-block:: shell

   git clone https://github.com/IMS-Bio2Core-Facility/BIC086-Sophie-Austin &&
   cd BIC086-Sophie-Austin &&
   snakemake --use-conda --cores 6

And ``snakemake`` will automatically leave it off.

Customising
~~~~~~~~~~~

You will almost certainly have a different gene list than me!
If its short,
you can specify it as a parameter like so:

.. code-block:: shell

   snakemake --use-conda --use-singularity --cores 6 --config gene_ids=["GENE1","GENE2"]

Alternatively, you can specify it as a list under the ``gene_ids`` parameter in the
configuration file located at ``configuration/snakemake.yaml``.

Unless you are querying a huge number of genes,
I find 6 cores to be sufficient to keep things moving.
If you want to use more or less,
just change the value passed to ``--cores``.
Each rule that implements concurrency will then use this many.
Remember, more isn't always faster.
On my laptop,
setting cores greater than the number of physical cores (not threads!)
in my machine gets me no improvement.
**YMMV**

Reproducibility
---------------

Reproducibility results are a cornerstone of the scientific process.
By running the pipeline with ``snakemake`` in a ``docker`` image using ``conda`` environments,
we ensure that no aspect of the pipeline is left to chance.
You will get our analysis,
as we ran it,
with the software versions,
as we used them.
To further aid in this effort,
`nox`_ and `pre-commit` are used,
which also ensures that development happens in reproducible environments.

Unfortunately,
any query to an unstable API is inherently not reproducible.
Thus,
changes in BioMart or GTEx could impact the results.
We recognise this as an inherent limitation,
and will do our best to keep abreast of API changes that impact the pipeline.

.. _nox: https://nox.thea.codes/en/stable/
.. _pre-commit: https://pre-commit.com/

Data
----

The pipeline requires no input data other than a list of gene names specified in
``configuration/snakemake.yaml``.

References
----------

It is surprisingly challenging to align RefSeq IDs and Ensembl IDs.
This is further complicated because GTEx uses Gencode26 under the hood.
As this is not the most up-to-date version,
it actually proved quite frustrating to find the desired version numbers for each gene.
To combat this,
this pipeline takes 3 different approaches in parallel:

#. Gencode v26 GTF annotations are downloaded from EBI,
   so the user only needs to supply gene names.
#. A query is made to BioMart to retrieve RefSeq IDs for each ENST returned by GTEx.
#. Data from `MANE`_ is added to help identify consensus transcripts.

.. _MANE: https://www.ncbi.nlm.nih.gov/refseq/MANE/

Contributing
------------

If you are interested in helping us improve the pipeline,
pleare see our guides on :ref:`contributing`.

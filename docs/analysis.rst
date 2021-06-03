.. _analysis:

Analysis
========

This analysis is performed inside the docker container
``condaforge/mambaforge:4.10.1-0``.

The entire ``Snakefile`` can be viewed at :ref:`Snakefile`.

Download Gencode v26 Reference
------------------------------

This step is handled entirely by the snakemake ``shell`` directive.
``wget`` is used to fetch the annotations,
which are then grepped to keep only those features marked ``gene``.
These features are then piped through a serious of ``cut`` and ``tr`` commands
to extract the relevant fields.

Download MANE Reference
-----------------------

.. automodule:: scripts.mane

Make GTEx Requests
------------------

.. automodule:: scripts.request

Make BioMart Requests
---------------------

.. automodule:: scripts.biomart

Process Final Data
------------------

.. automodule:: scripts.process

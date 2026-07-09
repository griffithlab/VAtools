VCF Expression Annotator
========================

The VCF Expression Annotator takes gene expression data and adds it to a VCF, allowing downstream tools to ask "how well is the gene containing this variant expressed?"

Usage
-----

.. program-output:: vcf-expression-annotator -h

Details
-----


This tool natively supports output files from `StringTie <https://github.com/gpertea/stringtie>`_, `Kallisto <https://kallisto.readthedocs.io/en/latest/>`_, or `Cufflinks <https://github.com/cole-trapnell-lab/cufflinks>`_, by specifying the appropriate format in the positional parameters: ``kallisto``, ``stringtie``, or ``cufflinks``.

In addition, the type of expression data, either ``gene`` or ``transcript``, needs to be specified. This will result in the expression value being written to the ``GX`` or ``TX`` field, respectively.

The input VCF needs to already be annotated with VEP, so that gene and transcript identifiers can be matched between the VCF and the expression file. 

When running in ``gene`` mode, Ensembl IDs - not gene names - are used. Depending on the expression software used, the Ensembl identifiers might contain version numbers. To add transcript and/or gene version numbers to your VEP annotations, use the ``--transcript_version`` and ``--gene-version`` when running VEP, respectively, as needed. 

You can also use the ``--ignore-ensembl-id-version`` flag of the VCF Expression Annotator to ignore the version of Ensembl gene and transcript IDs when finding the matching entry in your expression file.

### Custom Expression Data
VCF Expression Annotator can be used with other tools, so long as their expression output can be manipulated into a TSV containing two columns: Ensembl gene or transcript ID and expression values. This file also needs to contain a header line that is used to identify the contents of each column. These headers are specified via the ``--id-column`` and ``--expression-column`` parameters.

In order to use this option the ``custom`` value should be give in the file format parameter. Please note that when running in ``gene`` mode, the ID
column will need to contain Ensembl Gene IDs, not gene names.

By default the output VCF will be written to a ``.tx.vcf`` or ``.gx.vcf`` file next to
your input VCF file. You can set a different output file name using the
``--output-vcf`` parameter.

By default, running the VCF Expression Annotator against a VCF that already has a
``GX``/``TX`` FORMAT header raises an error. Use ``--force`` to annotate such a VCF
anyway - for example, to fill in a sample that wasn't previously annotated. Even with
``--force``, an error is raised if doing so would overwrite an existing non-blank
``GX``/``TX`` value for the target sample; use ``--overwrite`` (which requires
``--force``) to allow replacing those existing values.


Examples
--------

Download the example data used below:

.. code-block:: none

   curl -LO https://vatools.readthedocs.io/en/latest/_static/vatools-examples.tar.gz
   tar xzf vatools-examples.tar.gz
   cd vatools-examples

**1. Gene-level, adding Kallisto expression estimates from the gene-level TSV**

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.kallisto.genes.tsv kallisto gene \
     -o sample.gx.vcf

Adds a ``GX`` FORMAT field: ``ENSG00000184979|18.4008276190378``.

**2. Transcript-level, adding Kallisto expression estimates from the transcript-level TSV**

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.kallisto.transcripts.tsv kallisto transcript \
     -o sample.tx.vcf

Adds a ``TX`` FORMAT field: ``ENST00000215794|4``.

**3. StringTie, gene-level expression estimates from a TSV**

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.stringtie.genes.tsv stringtie gene \
     -o sample.gx.vcf

**4. StringTie, transcript-level expression estimates from a GTF**

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.stringtie.transcripts.gtf stringtie transcript \
     -o sample.tx.vcf

**5. Cufflinks, gene and transcript level**

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.cufflinks.genes.fpkm_tracking cufflinks gene \
     -o sample.gx.vcf

   vcf-expression-annotator sample.gx.vcf sample.cufflinks.isoforms.fpkm_tracking cufflinks transcript \
     -o sample.gx.tx.vcf

Chains gene- then transcript-level annotation onto the same VCF, so both FORMAT fields end up on one file.

**6.** Handle gene/tx ID suffixes with ``--ignore-ensembl-id-version``

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.kallisto.transcripts_versioned.tsv kallisto transcript \
     --ignore-ensembl-id-version -o sample.tx.vcf

The expression file's IDs carry a version suffix (the ``.1`` in ``ENST00000215794.1``) while the VCF's VEP annotation doesn't. Adding this flag allows matching. Without it, this command would produce a "missing expression entry" warning and no ``TX`` value.

**7. Using a custom format**

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.stringtie.genes.tsv custom gene \
     --id-column "Gene ID" --expression-column "TPM" -o sample.gx.vcf

Reuses the StringTie gene file to show that ``custom`` just needs ``--id-column``/``--expression-column`` pointed at whatever headers exist.

**8. Re-annotating with a different expression source (**``--force``/``--overwrite``**)**

.. code-block:: none

   vcf-expression-annotator sample.vcf sample.kallisto.genes.tsv kallisto gene \
     -o sample.gx.vcf

   vcf-expression-annotator sample.gx.vcf sample.stringtie.genes.tsv stringtie gene \
     --force --overwrite -o sample.gx.vcf

The first command annotates ``GX=ENSG00000184979|18.4008276190378`` from the Kallisto abundance file. Re-running against the same, already-annotated VCF with the StringTie TPM value for the same gene requires ``--force`` (to allow annotating a VCF that already has a ``GX`` header) and ``--overwrite`` on top of that (to replace the existing non-blank value) — the StringTie-derived ``GX=ENSG00000184979|2.629013`` overwrites the Kallisto-derived value.

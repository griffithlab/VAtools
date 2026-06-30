VCF Expression Annotator
========================

The VCF Expression Annotator takes gene expression data and adds it to a VCF, allowing downstream tools to ask "how well is the gene containing this variant expressed?"

The tool natively supports output files from [StringTie](https://github.com/gpertea/stringtie), [Kallisto](https://kallisto.readthedocs.io/en/latest/), or [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks), by specifying the appropriate format in the positional parameters: ``kallisto``, ``stringtie``, or ``cufflinks``.

In addition, the type of expression data, either ``gene`` or ``transcript``, needs to be specified. This will result in the expression value being written to the ``GX`` or ``TX`` field, respectively.

The input VCF needs to be annotated with VEP with gene and transcript information so that the VCF Expression Annotator can match a variant's Ensembl gene and transcript identifier in the VCF to the one in the expression file. 

When running in ``gene`` mode, Ensembl IDs - not gene names - are used. Depending on the expression software used, the transcript identifiers might contain version numbers. To add transcript version numbers to your VEP annotation, use the ``--transcript_version`` when running VEP. 

You can also use the ``--ignore-ensembl-id-version`` flag of the VCF Expression Annotator to ignore the version of Ensembl gene and transcript IDs when finding the matching entry in your expression file.

#### Custom Expression Data
VCF Expression Annotator can be used with other tools, so long as their expression output can be manipulated into a TSV containing two columns: Ensembl gene or transcript ID and expression values. This file also needs to contain a header line that is used to identify the contents of each column. These headers are specified via the ``--id-column`` and ``--expression-column`` parameters.

In order to use this option the ``custom`` value should be give in the file format parameter. Please note that when running in ``gene`` mode, the ID
column will need to contain Ensembl Gene IDs, not gene names.

By default the output VCF will be written to a ``.tx.vcf`` or ``.gx.vcf`` file next to
your input VCF file. You can set a different output file name using the
``--output-vcf`` parameter.

Usage
-----

.. program-output:: vcf-expression-annotator -h

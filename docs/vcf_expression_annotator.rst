VCF Expression Annotator
========================

The VCF Expression Annotator will take an output file from Cufflinks, Kallisto,
or StringTie and add the data from that file to your VCF. The expression file type is
specified using ``kallisto``, ``stringtie``, or ``cufflinks`` in the list of
positional parameters.

In addition, the type of expression data, either ``gene`` or ``transcript``, needs to
be specified. This will result in the expression value to be written to the
``GX`` or ``TX`` field, respectively.

The input VCF needs to be annotated with VEP with gene and transcript information so
that the VCF Expression Annotator can match a variant's gene and transcript
identifier in the VCF to the one in the expression file. Depending on the
expression software used, the transcript identifiers might contain version
numbers. To add transcript version numbers to your VEP annotation, use the
``--transcript_version`` when running VEP. You can also use the
``--ignore-ensembl-id-version`` flag of the VCF Expression Annotator to ignore
the version of Ensembl gene and transcript IDs when finding the matching entry in your expression
file.

The VCF Expression Annotator also accepts a custom tab-delimited (TSV) file input for the
expression file. This TSV file will need to contain one column with gene or
transcripts identifiers and one column with the expression values. This file
then needs to contain a header line that is used to
identify the contents of each column. This is done via the  ``--id-column``
and ``--expression-column`` parameters which need
to match the gene/transcript identifier and expression value column headers.
In order to use this option the expression file format option will need to be
set to ``custom``.

By default the output VCF will be written to a ``.tx.vcf`` or ``.gx.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.

Usage
-----

.. program-output:: vcf-expression-annotator -h

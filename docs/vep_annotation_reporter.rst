VEP Annotation Reporter
=======================

VEP annotations in a VCF can be difficult to read and interpret. We developed
the VEP Annotation Report to aid in converting VEP annotation fields to a
human-readable report.

The VEP Annotation Reporter will create a tab-delimited (TSV) file of
variants in a VCF and their VEP annotations. The VEP fields to add to the
output TSV are specified by listing one or more VEP annotation field identifiers as
positional parameters after the input VCF. The VEP fields of an
input VCF can be identified by inspecting the ``Description`` field of the
``CSQ`` ``INFO`` header in the VCF. Everything after ``Format:`` is a field
available in the VCF (delimited by ``|``).

If a variant is annotate with multiple transcript consequences by VEP then the
values for all transcript annotation will be returned as comma-separated
values. This is the default behavior unless VEP was run with
one of the ``--flag_pick`` options, all possible transcript consequences will be
reported by VEP but only one of these consequences will be picked by VEP as the
"best" consequence. This is denoted in the ``PICK`` field. If this field is
available, then the values for that transcript will be reported.

VEP annotations can also be added to an existing TSV with variant
information by using the ``--input-tsv`` option. In order to match
the variants in the TSV to the variants in the
VCF, the existing TSV file will need to contain columns with the headers
``CHROM``, ``POS``, ``REF``, and ``ALT`` where the values match the VCF
``CHROM``, ``POS``, ``REF``, and ``ALT`` values.

By default the output TSV will be written to a ``.tsv`` file next to
your input VCF file. You can set a different output file using the
``--output-tsv`` parameter.

Usage
-----

.. program-output:: vep-annotation-reporter -h

Example Command
---------------

.. code-block:: none

   vep-annotation-reporter input.vcf Consequence SYMBOL Feature -t input.tsv -o output.tsv

VEP Annotation Reporter
=======================

VEP annotations in a VCF are condensed into a CSQ field that is meant to be machine-readable, and can be difficult for humans to read and interpret. We developed the VEP Annotation Reporter to address this issue.

The VEP Annotation Reporter will create a tab-delimited (TSV) file of
variants in a VCF and their VEP annotations. The VEP fields to add to the
output TSV are specified by listing one or more VEP annotation field identifiers as
positional parameters after the input VCF. The VEP fields of an
input VCF can be identified by inspecting the ``Description`` field of the
``CSQ`` ``INFO`` header in the VCF. Everything after ``Format:`` is a field
available in the VCF (delimited by ``|``).

By default, VEP annotates each variant with it's consequences in every transcript it intersects. In this case, the values for all transcript annotation will be returned as comma-separated values. 

If VEP was run with one of the ``--flag_pick`` options, then it has labelled one of the consequences as the "best", and set the ``PICK`` field. If this field is
available, then only the values for that transcript will be reported. If no variant is set with the ``PICK`` field, then the values for all transcript consequences are reported.

VEP annotations can also be added to an existing TSV with variant
information by using the ``--input-tsv`` option. In order to match
the variants between the TSV and VCF, the existing TSV file will need to contain columns with the headers
``CHROM``, ``POS``, ``REF``, and ``ALT`` where the values exactly match the VCF
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

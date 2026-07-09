VEP Annotation Reporter
=======================

VEP annotations in a VCF are condensed into a CSQ field that is meant to be machine-readable, and can be difficult for humans to read and interpret. We developed the VEP Annotation Reporter to address this issue.

Usage
-----

.. program-output:: vep-annotation-reporter -h

Details
-----

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


Examples
--------

Download the example data used below:

.. code-block:: none

   curl -LO https://vatools.readthedocs.io/en/latest/_static/vatools-examples.tar.gz
   tar xzf vatools-examples.tar.gz
   cd vatools-examples

**1. Basic report**

.. code-block:: none

   vep-annotation-reporter sample.vcf Consequence SYMBOL Feature Amino_acids \
     -o sample.vep_report.tsv

Outputs one row: ``22 18644673 C T missense_variant USP18 ENST00000215794 A/V``.

**2. Multi-transcript variant, default (no PICK) behavior**

.. code-block:: none

   vep-annotation-reporter sample.multi_transcript.vcf Consequence SYMBOL Feature \
     -o sample.multi_report.tsv

Since VEP wasn't run with ``--flag_pick`` here, both transcripts' values are comma-separated in each field: ``intron_variant&nc_transcript_variant,non_coding_exon_variant&nc_transcript_variant``.

**3.** ``-p`` **preferred transcripts**

.. code-block:: none

   vep-annotation-reporter sample.multi_transcript.vcf Consequence SYMBOL Feature \
     -p sample.preferred_transcripts.list.tsv -o sample.preferred.tsv

Restricts the row to only ``ENST00000423297`` (``GRAMD4P2``) instead of both transcripts.

**4.** ``-t`` **merge into an existing TSV**

.. code-block:: none

   vep-annotation-reporter sample.vcf Consequence SYMBOL \
     -t sample.variant_report.tsv -o sample.merged.tsv

Appends ``Consequence``/``SYMBOL`` columns onto the existing ``sample.variant_report.tsv`` rows instead of writing a fresh CHROM/POS/REF/ALT TSV from scratch.

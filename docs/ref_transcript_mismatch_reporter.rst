Reference Transcript Mismatch Reporter
======================================

This tool can be used to identify variants where the reference genome build
doesn't match the Ensembl reference transcript used by VEP for variant consequence
annotations.

Usage
-----

.. program-output:: ref-transcript-mismatch-reporter -h

Details
-------

This situation is problematic because at these positions, the REF nucleotide(s) will differ
from the Ensembl transcript nucleotide(s) at the corresponding mutation
position. The end result is that any amino acid change predictions found in the ``Amino_acids``
field of the VEP CSQ annotation will be different from the translated Ensembl transcript
amino acids at that position.

This will lead to errors in some downstream tools (e.g. pVACseq) which rely on accurate peptide sequences in order to make predictions about neoantigen MHC binding. When mismatches occur, these two
fields will be in disagreement and pVACseq
cannot make predictions on such variants.

If there are only minor differences between the reference used and the Ensembl transcripts, only a small number of variants may be affected. If a large number of variants are flagged, then it would be wise to step back and consider whether there is a fundamental lack of compatibility (for example, alignments to GRCh37 that use a GRCh38 VEP cache).

**Inputs**

The input VCF needs to be annotated by VEP and requires annotation with the
``Wildtype`` VEP plugin available as part of `pVACtools <https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html#installing-vep>`_.

**Outputs**

This tool will report on the number of variants and transcripts in a VCF that
are affected by this issue and output this information to stdout. It will
write a ``.mismatch.tsv`` file next to the VCF that provides further details
on the problematic variants.

**Filtering**

This tool also allows the user to either soft-filter or hard-filter the VCF
using the ``--filter [soft|hard]`` parameter. Soft-filtering will tag the
problematic variants with a custom VCF FILTER ``CSQ_MISMATCH`` while hard-filtering
will produce a new VCF that has these variants removed. When using a filter,
the output VCF will be written to a ``filtered.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.


Examples
--------

Download the example data used below:

.. code-block:: none

   curl -LO https://vatools.readthedocs.io/en/latest/_static/vatools-examples.tar.gz
   tar xzf vatools-examples.tar.gz
   cd vatools-examples

These examples use ``sample.mismatch.vcf``, which was deliberately constructed so the REF base disagrees with the Ensembl transcript sequence VEP used.

**1. Report only**

.. code-block:: none

   ref-transcript-mismatch-reporter sample.mismatch.vcf

Writes ``sample.mismatches.tsv`` next to the input and prints summary counts (variants processed, percentage flagged) to stdout.

**2. Soft-filter**

.. code-block:: none

   ref-transcript-mismatch-reporter sample.mismatch.vcf --filter soft \
     -o sample.soft_filtered.vcf

Tags mismatched records with ``FILTER=CSQ_MISMATCH`` instead of removing them.

**3. Hard-filter**

.. code-block:: none

   ref-transcript-mismatch-reporter sample.mismatch.vcf --filter hard \
     -o sample.hard_filtered.vcf

Removes mismatched records from the output VCF entirely.

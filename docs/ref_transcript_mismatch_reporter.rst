Reference Transcript Mismatch Reporter
======================================

This tool can be used to identify variants where the reference genome build
doesn't match the Ensembl reference transcript used by VEP for variant consequence
annotations. This is problematic because at these positions the REF nucleotide(s) will differ
from the Ensembl transcript nucleotide(s) at the corresponding mutation
position. The end result is that any amino acid change predictions found in the ``Amino_acids``
field of the VEP CSQ annotation will be different from the translated Ensembl transcript
amino acids at that position.

This will lead to errors in some downstream tools (e.g. pVACseq) which rely on accurate peptide sequences in order to make predictions about neoantigen MHC binding. When mismatches occur, these two
fields will be in disagreement and pVACseq
cannot make predictions on such variants.

If there are only minor differences between the reference used and the Ensembl transcripts, only a small number of variants may be affected. If a large number of variants are flagged, then it would be wise to step back and consider whether there is a fundamental lack of compatibility (for example, alignments to GRCh37 that use a GRCh38 VEP cache).

The input VCF needs to be annotated by VEP and requires annotation with the
``Wildtype`` VEP plugin available as part of `pVACtools <https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html#installing-vep>`_.

This tool will report on the number of variants and transcripts in a VCF that
are affected by this issue and output this information to stdout. It will
write a ``.mismatch.tsv`` file next to the VCF that provides further details
on the problematic variants.

This tool also allows the user to either soft-filter or hard-filter the VCF
using the ``--filter [soft|hard]`` parameter. Soft-filtering will tag the
problematic variants with a custom VCF FILTER ``CSQ_MISMATCH`` while hard-filtering
will produce a new VCF that has these variants removed. When using a filter,
the output VCF will be written to a ``filtered.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.

Usage
-----

.. program-output:: ref-transcript-mismatch-reporter -h

Example Command
---------------

.. code-block:: none

   ref-transcript-mismatch-reporter input.vcf --filter soft
